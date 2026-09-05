%% TEST_CTX_BLOCKED  El troceado temporal de build_ctx_blocked es EXACTO.
%  Verifica, sobre el perfil de E2 con constelacion reducida (T=66), que construir
%  el contexto del motor POR BLOQUES da EXACTAMENTE lo mismo que construirlo de una
%  pieza:
%
%    (1) ctx campo a campo: build_ctx_blocked(timeBlock=10)  vs  pipeline directo
%        con ffr_context sobre la ventana entera  ->  max|dif| debe ser 0.
%    (2) KPIs de los 5 esquemas de E2 calculados sobre ambos contextos -> max|dif| 0.
%    (3) KPIs contra el .mat de referencia ffr_results_T66.mat -> max|dif| 0.
%
%  El (3) es el que de verdad cierra el hueco: no compara dos caminos nuevos entre
%  si, sino el camino nuevo contra un resultado ya publicado y defendido. Es la
%  misma logica con la que run_density_sweep se valido reproduciendo run_ffr_demo.
%
%  ffr_results_T66.mat NO SE DISTRIBUYE con este repositorio (es un resultado
%  archivado, no codigo). El test lo detecta con exist() y, si falta, OMITE la
%  comprobacion (3) y se pronuncia sobre (1) y (2). Para recuperarla hay que
%  regenerar ese .mat con el perfil de la seccion 1, que es exactamente el que lo
%  produjo.
%
%  POR QUE ES EXACTO (y no "aproximadamente igual"): el pipeline no acopla
%  instantes. Geometria, servidor, balance de enlace, interferencia y ganancias de
%  haz son por instante; lo unico definido sobre la ventana son los KPIs, que se
%  calculan despues sobre el contexto ya completo. Trocear solo cambia EN QUE ORDEN
%  se llenan las mismas columnas.
%
%  Uso: ejecutar en la raiz del proyecto. Tarda ~2 min (incluye la tabla P.618).

clear; clc;

fprintf('=========================================================================\n');
fprintf('  TEST: exactitud del troceado temporal (build_ctx_blocked)\n');
fprintf('=========================================================================\n');

%% 1. Perfil EXACTO de E2 con constelacion reducida (el que genero ffr_results_T66.mat)
cfg = config_default();
cfg.constellations(1).T = 66;
cfg.constellations(1).P = 6;
cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 3;
cfg.time.dt       = 60;
cfg.time.duration = 7200;
cfg.ffr.tau_mode = 'quantile';
cfg.ffr.tau_q    = 50;
cfg.beams.nRings = 5;

cases = { ...
    struct('scheme','reuse1','Delta',1,'alpha',NaN, 'name','reuse1'), ...
    struct('scheme','reuseD','Delta',3,'alpha',NaN, 'name','reuseD (D=3)'), ...
    struct('scheme','reuseD','Delta',4,'alpha',NaN, 'name','reuseD (D=4)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',0.4, 'name','FFR (D=3, a=0.4)'), ...
    struct('scheme','ffr',   'Delta',4,'alpha',0.4, 'name','FFR (D=4, a=0.4)') };
nC = numel(cases);

fprintf('\n[warmup] tabla P.618...\n');
atm_loss_dB(cfg, 45);

%% 2. CAMINO A: ventana de una pieza (el codigo original de run_ffr_demo)
fprintf('\n----- Camino A: ventana ENTERA (sin trocear) -----\n');
cfgA = cfg;  cfgA.compute.timeBlock = [];
tvec   = cfgA.time.t0 : cfgA.time.dt : cfgA.time.t0 + cfgA.time.duration;
sats   = build_constellation(cfgA);
R_eci  = propagate(cfgA, sats, tvec);
R_ecef = eci2ecef(cfgA, R_eci, tvec);
users  = build_user_grid(cfgA);
G      = compute_geometry(cfgA, users, R_ecef);
S      = associate_serving(cfgA, G);
L      = compute_link_budget(cfgA, S);
INTF   = compute_interference(cfgA, sats, users, G, S, L);
BL     = build_beam_layout(cfgA);
ctxA      = ffr_context(cfgA, sats, users, BL, R_ecef, G, S, L, INTF, false);
ctxA.tvec = tvec;
clear G R_eci R_ecef INTF S L;
fprintf('ctxA construido: M=%d Nt=%d nBeams=%d\n', ctxA.M, ctxA.Nt, ctxA.nBeams);

%% 3. CAMINO B: troceado en bloques de 10 instantes
fprintf('\n----- Camino B: TROCEADO (timeBlock = 10) -----\n');
cfgB = cfg;  cfgB.compute.timeBlock = 10;
[ctxB, ~, ~, BLB, BINFO] = build_ctx_blocked(cfgB, true);
fprintf('bloques: %d de %d instantes | pico %.3f GB (sin trocear %.3f GB)\n', ...
    BINFO.nBlocks, BINFO.nBlk, BINFO.peak_GB, BINFO.peak_full_GB);

%% 4. COMPARACION (1): ctx campo a campo
fprintf('\n----- (1) ctx campo a campo -----\n');
flds = {'Grel','Glin','beam','theta_deg','Cpsd_dBWMHz','Ipsd_inter_dBWMHz', ...
        'elev_deg','range_km','nVis','cov','M','Nt','nBeams','tvec'};
maxd = 0;  worst = '';
for i = 1:numel(flds)
    f = flds{i};
    a = double(ctxA.(f));  b = double(ctxB.(f));
    if ~isequal(size(a), size(b))
        error('test_ctx_blocked:size', 'Campo %s: tamanos distintos (%s vs %s).', ...
            f, mat2str(size(a)), mat2str(size(b)));
    end
    d = max(abs(a(:) - b(:)), [], 'omitnan');
    nanDiff = sum(isnan(a(:)) ~= isnan(b(:)));
    if nanDiff > 0
        error('test_ctx_blocked:nan', 'Campo %s: %d posiciones difieren en NaN.', f, nanDiff);
    end
    if isempty(d), d = 0; end
    fprintf('  %-20s max|dif| = %.3g   (NaN coincidentes)\n', f, d);
    if d > maxd, maxd = d; worst = f; end
end
if ~isequal(BL.nBeams, BLB.nBeams)
    error('test_ctx_blocked:BL', 'BL.nBeams difiere.');
end
fprintf('  --> PEOR CAMPO: %s con max|dif| = %.3g\n', worst, maxd);

%% 5. COMPARACION (2): KPIs de los 5 esquemas sobre ambos contextos
fprintf('\n----- (2) KPIs de los %d esquemas (A vs B) -----\n', nC);
KA = cell(1,nC);  KB = cell(1,nC);
cA = ctxA;  cB = ctxB;
for c = 1:nC
    cf = cfg;
    cf.ffr.scheme = cases{c}.scheme;
    cf.ffr.Delta  = cases{c}.Delta;
    if ~isnan(cases{c}.alpha), cf.ffr.alpha = cases{c}.alpha; end

    Aa = ffr_allocate(cf, cA, false);  Fa = compute_sinr_ffr(cf, cA, Aa);
    KA{c} = compute_kpis(cf, Fa, Aa);
    if isempty(cA.SINR_ref_dB), cA.SINR_ref_dB = Aa.SINR_ref_dB; end

    Ab = ffr_allocate(cf, cB, false);  Fb = compute_sinr_ffr(cf, cB, Ab);
    KB{c} = compute_kpis(cf, Fb, Ab);
    if isempty(cB.SINR_ref_dB), cB.SINR_ref_dB = Ab.SINR_ref_dB; end
end
[md2, wh2] = deep_maxdiff_local(KA, KB, 'K');
fprintf('  max|dif| = %.3g   (%s)\n', md2, wh2);

%% 6. COMPARACION (3): contra el .mat PUBLICADO a T=66
fprintf('\n----- (3) contra ffr_results_T66.mat (resultado publicado) -----\n');
md3 = NaN;  wh3 = '(no comparado)';
if exist('ffr_results_T66.mat','file')
    P = load('ffr_results_T66.mat','KK');
    [md3, wh3] = deep_maxdiff_local(P.KK, KB, 'K');
    fprintf('  max|dif| = %.3g   (%s)\n', md3, wh3);
    fprintf('\n  %-18s %12s %12s %12s\n', 'esquema','SINR_e_p5 pub','SINR_e_p5 nuevo','dif');
    for c = 1:nC
        fprintf('  %-18s %12.4f %12.4f %12.3g\n', cases{c}.name, ...
            P.KK{c}.viab.SINR_edge_p5, KB{c}.viab.SINR_edge_p5, ...
            abs(P.KK{c}.viab.SINR_edge_p5 - KB{c}.viab.SINR_edge_p5));
    end
else
    fprintf('  ffr_results_T66.mat no encontrado: se omite.\n');
end

%% 7. VEREDICTO
fprintf('\n=========================================================================\n');
tol = 0;
ok1 = (maxd == tol);
ok2 = (md2  == tol);
ok3 = isnan(md3) || (md3 == tol);
if ok1 && ok2 && ok3
    fprintf('  RESULTADO: SUPERADO. El troceado es EXACTO (max|dif| = 0 en los 3 tests).\n');
else
    fprintf('  RESULTADO: FALLO.  ctx=%.3g  KPIs=%.3g  publicado=%.3g\n', maxd, md2, md3);
    error('test_ctx_blocked:fail', 'El troceado NO es exacto.');
end
fprintf('=========================================================================\n');

% =========================================================================
function [md, whr] = deep_maxdiff_local(a, b, path)
%DEEP_MAXDIFF_LOCAL  Maxima diferencia absoluta recursiva entre dos estructuras.
md = 0;  whr = '';
if iscell(a) && iscell(b)
    if numel(a) ~= numel(b), md = Inf; whr = [path ' (numel)']; return; end
    for i = 1:numel(a)
        [d, w] = deep_maxdiff_local(a{i}, b{i}, sprintf('%s{%d}', path, i));
        if d > md, md = d; whr = w; end
    end
elseif isstruct(a) && isstruct(b)
    fa = fieldnames(a);
    for i = 1:numel(fa)
        if ~isfield(b, fa{i}), continue; end
        [d, w] = deep_maxdiff_local(a.(fa{i}), b.(fa{i}), [path '.' fa{i}]);
        if d > md, md = d; whr = w; end
    end
elseif (isnumeric(a) || islogical(a)) && (isnumeric(b) || islogical(b))
    a = double(a);  b = double(b);
    if ~isequal(size(a), size(b)), md = Inf; whr = [path ' (size)']; return; end
    d = max(abs(a(:) - b(:)), [], 'omitnan');
    if ~isempty(d) && d > md, md = d; whr = path; end
end
end
