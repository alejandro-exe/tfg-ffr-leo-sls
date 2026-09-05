function [ctx, sats, users, BL, INFO] = build_ctx_blocked(cfg, verbose)
%BUILD_CTX_BLOCKED  Contexto del motor FFR construido por BLOQUES TEMPORALES.
%   [ctx, sats, users, BL, INFO] = build_ctx_blocked(cfg)
%   [ctx, sats, users, BL, INFO] = build_ctx_blocked(cfg, verbose)
%
%   HABILITADOR DE MEMORIA de los runners de un solo punto (E2 = run_ffr_demo y
%   E5 = run_e5_adaptive) a densidad ORBITAL REAL. Recorre la ventana temporal en
%   bloques de cfg.compute.timeBlock instantes ejecutando el pipeline dependiente
%   del tiempo
%
%       propagate -> eci2ecef -> compute_geometry -> associate_serving
%                 -> compute_link_budget -> compute_interference -> ffr_context
%
%   y CONCATENA los bloques en un unico ctx identico al que devolveria ffr_context
%   sobre la ventana entera. Los runners no cambian aguas abajo: siguen recibiendo
%   el ctx completo y todo su analisis (barrido de alpha, desglose por elevacion,
%   p5split, figuras) funciona sin tocar una linea.
%
%   POR QUE HACE FALTA. El pico de memoria del pipeline es [M x N x Ntb] dentro de
%   compute_interference (~74 B/elemento, ver mem_peak_model_GB). Con el perfil de
%   E2/E5 (M=553 usuarios, Nt=121 instantes, 91 haces) y la constelacion REAL
%   T=1584 son 7.53 GB de una pieza, MAS que la RAM total de la maquina de
%   referencia (7.45 GB): sin trocear, los experimentos de cabecera simplemente no
%   se pueden ejecutar a densidad real. Con timeBlock = 10 el pico lo fija el
%   BLOQUE y baja a ~0.62 GB.
%
%   POR QUE ES EXACTO (y no una aproximacion). El pipeline NO acopla instantes: la
%   geometria, la asociacion de servidor, el balance de enlace, la interferencia y
%   las ganancias de haz se calculan por instante. Lo unico definido sobre la
%   ventana son los KPIs (percentiles y medias temporales), que se calculan DESPUES
%   sobre el ctx ya completo. Es el mismo argumento -- y el mismo troceado -- que
%   run_one_density usa en los barridos, verificado en test_timeblock_invariance
%   (max|dif| = 0) y aqui de nuevo contra el .mat publicado de E2 a T=66.
%
%   FUENTE UNICA. E2 y E5 comparten esta funcion en lugar de llevar cada uno su
%   copia del bucle: es la politica del proyecto (merge_cfg, mem_peak_model_GB,
%   grid_outside_cluster) y la que evita el modo de fallo que ya causo los bugs de
%   atm_key_of y associate_serving -- la misma regla escrita dos veces que acaba
%   divergiendo sin que nadie lo note.
%
%   Entradas:
%     cfg     - configuracion completa (usa cfg.compute.timeBlock; [] / Inf = sin
%               trocear, un solo bloque, camino de codigo historico).
%     verbose - true (defecto): imprime el diagnostico de ffr_context UNA vez sobre
%               el ctx COMPLETO (no una vez por bloque).
%
%   Salidas:
%     ctx   - contexto completo [M x Nt] / [M x nBeams x Nt], con ctx.tvec relleno.
%     sats, users, BL - construidos aqui (no dependen del instante).
%     INFO  - .nBlk .nBlocks .peak_GB (bloque) .peak_full_GB (sin trocear) .M .N
%             .Nt .nBeams .elapsed_s
%
%   Aportacion propia (infraestructura; NO contiene fisica ni altera ningun
%   resultado: reordena EN EL TIEMPO las mismas llamadas al mismo motor).

if nargin < 2 || isempty(verbose), verbose = true; end

tB = tic;

%% 1. Piezas que NO dependen del instante
tvec  = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
Nt    = numel(tvec);
sats  = build_constellation(cfg);
users = build_user_grid(cfg);
BL    = build_beam_layout(cfg);
M     = users.M;   N = sats.N;   nB = BL.nBeams;

% Satelites ELEGIBLES como servidores (cfg.geom.serving_cid = [] -> todos, que es
% el defecto historico). Mismo criterio que run_one_density.
servMask = [];
if isfield(cfg,'geom') && isfield(cfg.geom,'serving_cid') && ~isempty(cfg.geom.serving_cid)
    servMask = ismember(sats.cid, cfg.geom.serving_cid);
    if ~any(servMask)
        error('build_ctx_blocked:servingCid', ...
            'cfg.geom.serving_cid = [%s] no selecciona ningun satelite.', ...
            num2str(cfg.geom.serving_cid));
    end
end

%% 2. Tamano de bloque (MISMA regla que run_one_density: lectura defensiva)
tb = Inf;
if isfield(cfg,'compute') && isfield(cfg.compute,'timeBlock') && ~isempty(cfg.compute.timeBlock)
    tb = cfg.compute.timeBlock;
end
if ~isfinite(tb) || tb <= 0 || tb >= Nt
    nBlk = Nt;
else
    nBlk = max(1, floor(tb));
end
nBlocks = ceil(Nt / nBlk);

peakBlk  = mem_peak_model_GB(M, N, nBlk, nB);
peakFull = mem_peak_model_GB(M, N, Nt,   nB);

if verbose
    fprintf(['[troceado] %d instantes en %d bloque(s) de %d | pico modelado %.2f GB ' ...
             '(sin trocear: %.2f GB, x%.2f)\n'], Nt, nBlocks, nBlk, peakBlk, peakFull, ...
             peakFull/max(peakBlk,eps));
end

%% 3. Acumuladores del ctx COMPLETO
%  Ocupan M*nBeams*Nt (decenas de MB), tres ordenes por debajo del M*N*Ntb que se
%  quiere acotar: acumular el ctx entero NO compromete el objetivo del troceado.
ctx        = struct();
ctx.BL     = BL;
ctx.M      = M;
ctx.Nt     = Nt;
ctx.nBeams = nB;
ctx.Grel   = nan(M, nB, Nt);
ctx.Glin   = nan(M, nB, Nt);
ctx.beam      = nan(M, Nt);
ctx.theta_deg = nan(M, Nt);
ctx.Cpsd_dBWMHz       = nan(M, Nt);
ctx.Ipsd_inter_dBWMHz = nan(M, Nt);
ctx.elev_deg = nan(M, Nt);
ctx.range_km = nan(M, Nt);
ctx.nVis     = nan(M, Nt);
ctx.cov      = false(M, Nt);

%% 4. Bucle de bloques
for b = 1:nBlocks
    kIdx = ((b-1)*nBlk + 1) : min(b*nBlk, Nt);
    tv   = tvec(kIdx);

    R_eci  = propagate(cfg, sats, tv);
    R_ecef = eci2ecef(cfg, R_eci, tv);
    % sats.minElev = mascara POR CONSTELACION; con una sola constelacion y
    % minElev = [] equivale al escalar global (bit-identico, ver compute_geometry).
    G      = compute_geometry(cfg, users, R_ecef, sats.minElev);
    S      = associate_serving(cfg, G, servMask);
    L      = compute_link_budget(cfg, S);
    INTF   = compute_interference(cfg, sats, users, G, S, L);

    % verbose = false SIEMPRE: el diagnostico de ffr_context es sobre la ventana y
    % se emite una sola vez al final, sobre el ctx ya completo (mas abajo).
    cb = ffr_context(cfg, sats, users, BL, R_ecef, G, S, L, INTF, false);

    % Liberar AQUI es lo que ACOTA el pico en vez de solo repartirlo: si G siguiera
    % viva, el bloque siguiente sumaria su geometria a la anterior.
    clear G R_eci R_ecef INTF S L;

    ctx.Grel(:,:,kIdx) = cb.Grel;
    ctx.Glin(:,:,kIdx) = cb.Glin;
    ctx.beam(:,kIdx)      = cb.beam;
    ctx.theta_deg(:,kIdx) = cb.theta_deg;
    ctx.Cpsd_dBWMHz(:,kIdx)       = cb.Cpsd_dBWMHz;
    ctx.Ipsd_inter_dBWMHz(:,kIdx) = cb.Ipsd_inter_dBWMHz;
    ctx.elev_deg(:,kIdx) = cb.elev_deg;
    ctx.range_km(:,kIdx) = cb.range_km;
    ctx.nVis(:,kIdx)     = cb.nVis;
    ctx.cov(:,kIdx)      = cb.cov;
    clear cb;
end

ctx.tvec        = tvec;
ctx.SINR_ref_dB = [];        % lo rellena ffr_allocate / el runner (igual que ffr_context)

%% 5. Diagnostico sobre el ctx COMPLETO (mismo texto que ffr_context, una sola vez)
if any(ctx.cov(:))
    thv  = ctx.theta_deg(ctx.cov);
    latt = 'hex';
    if isfield(BL,'lattice') && ~isempty(BL.lattice), latt = BL.lattice; end
    [isOut, pctCirc, thLim] = grid_outside_cluster(max(thv), BL.spacing_deg, latt);
    if strcmp(latt,'rect'), circLbl = 's*sqrt(2)/2'; else, circLbl = 's/sqrt(3)'; end

    if verbose
        fprintf('[ffr_context] M=%d usuarios | Nt=%d instantes | %d haces\n', M, Nt, nB);
        fprintf(['[ffr_context] offset angular al haz servidor: media %.3f deg, max %.3f deg ' ...
                 '(s=%.3f deg | circunradio %s=%.3f deg -> %.1f%%)\n'], ...
            mean(thv), max(thv), BL.spacing_deg, circLbl, thLim, pctCirc);
        fprintf('[ffr_context] elevacion del servidor: media %.1f deg, min %.1f deg, max %.1f deg\n', ...
            mean(ctx.elev_deg(ctx.cov)), min(ctx.elev_deg(ctx.cov)), max(ctx.elev_deg(ctx.cov)));
    end

    if isOut
        warning('ffr_context:offset', ...
            ['Hay usuarios FUERA de su celda (max offset %.3f deg = %.1f%% del ' ...
             'circunradio de celda %s = %.3f deg, con s = %.3f deg): la rejilla de ' ...
             'usuarios sobresale del cluster de haces.'], ...
            max(thv), pctCirc, latt, thLim, BL.spacing_deg);
    end
end

INFO = struct('nBlk', nBlk, 'nBlocks', nBlocks, 'peak_GB', peakBlk, ...
              'peak_full_GB', peakFull, 'M', M, 'N', N, 'Nt', Nt, ...
              'nBeams', nB, 'elapsed_s', toc(tB));
end
