function OUT = run_one_density(cfg_base, pt)
%RUN_ONE_DENSITY  UN punto del barrido de saturacion, de extremo a extremo.
%   OUT = run_one_density(cfg_base, pt)
%
%   FUNCION PURA DE PUNTO (preparacion de E3 / Monte Carlo). Ejecuta el pipeline
%   completo para una configuracion de densidad y devuelve SOLO structs de
%   resultados. No lee ni escribe estado global, no dibuja, no guarda ficheros y
%   no modifica cfg_base: trabaja sobre una COPIA local. Por eso dos puntos del
%   barrido son estrictamente independientes y el bucle que los recorre puede ser
%   `for` o `parfor` indistintamente (ver run_sweep_points).
%
%   PIPELINE (el mismo que run_ffr_demo; aqui NO hay fisica nueva):
%     build_constellation -> propagate -> eci2ecef -> build_user_grid
%       -> compute_geometry -> associate_serving -> compute_link_budget
%       -> compute_interference -> build_beam_layout -> ffr_context
%       -> [por esquema] ffr_allocate -> compute_sinr_ffr -> compute_kpis
%
%   ---------------------------------------------------------------------------
%   TROCEADO TEMPORAL (cfg.compute.timeBlock) -- control del PICO de RAM
%   ---------------------------------------------------------------------------
%   El pico de memoria de un punto NO son los tensores de compute_geometry: se
%   alcanza dentro de compute_interference, en atm_loss_dB(cfg, G.el), donde sobre
%   un array [M x N x Nt] coexisten G.el/az/range (24 B/elem) + G.vis (1) +
%   FSPL_all (8) + Latm_clear/rain (16) + valid (1) + temporales de interp1 (~24),
%   es decir ~74 B por elemento M*N*Nt (modelo de la FASE A; medido 61 B, luego es
%   cota superior). Con 317 usuarios, T=1584 y Nt=61 eso son ~2.2 GB por punto, y a
%   T=4000 son ~5.4 GB: mas de lo que cabe en un portatil, y multiplicado por cada
%   worker de parfor.
%
%   Con cfg.compute.timeBlock = n la ventana se procesa en BLOQUES de n instantes
%   y el pico pasa a fijarlo n, NO Nt:  pico ~ 74*M*N*n.  Es una optimizacion de
%   RENDIMIENTO/MEMORIA, no un cambio de modelo:
%     - el pipeline NO acopla instantes distintos (geometria, radioenlace e
%       interferencia son por instante; ffr_policy fija tau por instante en modo
%       cuantil y alpha por instante; ffr_allocate mide la carga n0 por instante);
%     - y aqui NO se combinan KPIs parciales, que es donde seria facil equivocarse
%       (el percentil 5 global NO es la media de los p5 de cada bloque). Se acumulan
%       las MATRICES DE MUESTRA [M x Nt] -- SINR, B_user, isCenter, C/N,
%       penalizacion -- y se llama a compute_kpis UNA sola vez al final sobre la
%       matriz completa, que es literalmente el mismo array que recibiria sin
%       trocear. Percentiles, CDF, Jain, R_agg y el veredicto de viabilidad salen
%       por tanto EXACTOS por construccion, no "aproximadamente iguales".
%   Los acumuladores ocupan M*Nt (cientos de kB), despreciable frente a M*N*Nt.
%   Verificado en test_timeblock_invariance: troceado vs sin trocear, max|dif| = 0.
%
%   cfg.compute.timeBlock = [] o Inf  ->  un solo bloque = comportamiento de antes.
%
%   PORQUE VARIOS ESQUEMAS POR PUNTO: la geometria (M*N*Nt) y ffr_context son el
%   ~90% del coste y NO dependen del esquema. Evaluar aqui los cuatro esquemas de
%   E3 (reuse1 / reuseD / ffr-estatica / ffr-adaptativa) sobre el MISMO ctx es a
%   la vez mas rapido y metodologicamente obligatorio: comparten poblacion de
%   usuarios, instantes y clasificacion centro/borde de referencia.
%
%   Entradas:
%     cfg_base - configuracion base (p.ej. config_default con el perfil de estudio).
%                NO se modifica: se copia por valor.
%     pt       - struct del PUNTO:
%       .T, .P            (opc.) sobrescriben cfg.constellations(1) [densidad]
%       .constellations   (opc.) reemplaza la LISTA entera de constelaciones
%                                (caso INTER-constelacion de la saturacion)
%       .cfg_over         (opc.) struct ANIDADO que se funde sobre cfg (merge
%                                recursivo: solo pisa las hojas que trae). Es el
%                                gancho para barrer CUALQUIER eje distinto de la
%                                densidad orbital sin tocar esta funcion; E3 lo usa
%                                para el eje de DENSIDAD DE HACES
%                                (.radio.beamwidth3dB_deg / .beams.nRings / la
%                                EIRP acoplada). run_sweep_points YA lo asumia al
%                                calcular la clave de la cache P.618 (atm_key_of),
%                                asi que esto cierra un hueco, no abre uno nuevo.
%                                Se aplica ANTES que .T/.P/.constellations.
%       .cases            (opc.) cell de structs de esquema, cada uno con campos
%                                {scheme, Delta, alpha, adaptive, name} (los que
%                                falten se dejan como esten en cfg_base.ffr)
%       .scheme/.Delta/.alpha/.adaptive   alternativa a .cases para UN solo esquema
%       .elMin_valid      (opc.) deg: ademas de la ventana completa, KPIs sobre la
%                                submuestra de GEOMETRIA VALIDA (elev >= este valor)
%       .idx              (opc.) indice del punto -> substream del RNG (ver abajo)
%       .seed             (opc.) semilla base del punto (defecto 42)
%       .keep_cdf         (opc.) true = conservar las CDF en los KPIs (defecto false:
%                                son vectores de decenas de miles de puntos y en un
%                                barrido inflan la transferencia worker -> cliente)
%       .label            (opc.) etiqueta de texto del punto
%       .verbose          (opc.) true = log del punto (defecto false)
%
%   Salida (struct OUT, TODO son structs/escalares pequenos):
%     .label, .T, .P, .N, .M, .Nt, .nBeams, .Delta_list
%     .cases      cell con la definicion de cada esquema evaluado
%     .K          cell de KPIs (compute_kpis) en la VENTANA COMPLETA
%     .KV         cell de KPIs sobre la submuestra de geometria valida ([] si no)
%     .diag       diagnostico del punto (cobertura, elevaciones, criterio de H2)
%     .mem        troceado temporal aplicado y coste de memoria (informativo):
%                 .timeBlock/.nBlocks, .peak_model_GB (con troceado),
%                 .peak_full_GB (de una pieza), .samples_MB (acumuladores)
%     .runinfo    metadatos NO deterministas (tiempo de calculo, worker): se
%                 excluyen EXPRESAMENTE de cualquier comparacion de invariancia.
%
%   ---------------------------------------------------------------------------
%   REPRODUCIBILIDAD (obligatoria para que parfor == for)
%   ---------------------------------------------------------------------------
%   El pipeline de la FFR es DETERMINISTA: no hay ninguna llamada a rand/randn en
%   la cadena (geometria kepleriana, patrones analiticos, P.618 tabulada,
%   percentiles). El unico rand del proyecto vive en compute_beam_cir /
%   run_calibration_e0 (E0, usuarios sinteticos), que NO forma parte de este
%   pipeline y ademas fija rng(0) internamente.
%
%   Aun asi se fija AQUI un stream de numeros aleatorios PROPIO DEL PUNTO
%   (RandStream con Substream = pt.idx), no global y no dependiente del orden de
%   ejecucion de los workers. Motivo: en cuanto el Monte Carlo de E3 introduzca
%   sorteos (fases orbitales, posiciones de usuario, desvanecimientos), el
%   resultado del punto i debe seguir siendo el MISMO tanto si se calcula el
%   primero como el ultimo, en el cliente o en cualquier worker. Con semilla
%   global (rng(s) en el runner) eso NO se cumple bajo parfor.
%
%   ---------------------------------------------------------------------------
%   CACHE PERSISTENT DE atm_loss_dB BAJO PARFOR
%   ---------------------------------------------------------------------------
%   atm_loss_dB cachea la tabla P.618 en una variable `persistent`, cuya clave es
%       freq_GHz | ground.point(1) | ground.point(2) | geom.minElev | p618_availability
%   Cada worker de parfor es un proceso MATLAB independiente con su propio espacio
%   de variables persistent: construira la tabla UNA VEZ (~90 s) y la reutilizara
%   en todos los puntos que le toquen. Es el comportamiento correcto y deseado.
%   NINGUNO de los parametros de la clave depende de la densidad de satelites, asi
%   que un barrido de saturacion completo usa UNA sola tabla por worker. Si un
%   futuro barrido moviese la frecuencia, el emplazamiento o la disponibilidad,
%   la clave cambiaria y la tabla se reconstruiria sola: la cache NUNCA devuelve
%   valores de otra configuracion (verificado en test_atm_cache.m).
%
%   Aportacion propia (infraestructura de experimentos; no altera la fisica).

%% 0. Copia LOCAL de la configuracion (nada compartido entre puntos)
cfg = cfg_base;

% Overrides genericos de cfg (merge recursivo). Van PRIMERO para que .T/.P y
% .constellations sigan mandando sobre la densidad orbital si se dan ambos.
if isfield(pt,'cfg_over') && isstruct(pt.cfg_over) && ~isempty(fieldnames(pt.cfg_over))
    cfg = merge_cfg(cfg, pt.cfg_over);
end

if isfield(pt,'constellations') && ~isempty(pt.constellations)
    cfg.constellations = pt.constellations;
end
if isfield(pt,'T') && ~isempty(pt.T), cfg.constellations(1).T = pt.T; end
if isfield(pt,'P') && ~isempty(pt.P), cfg.constellations(1).P = pt.P; end

% Coherencia de la cfg del PUNTO, justo despues de aplicar los overrides y ANTES de
% simular nada: es donde pt.cfg_over puede haber pisado un maestro (B_MHz, GT_dBK,
% Grx_max_dBi) dejando su derivado obsoleto. Fallar aqui cuesta milisegundos; fallar
% al final de un barrido cuesta minutos y da numeros plausibles pero equivocados.
% (auditoria, hallazgo G3; el chequeo se repite en compute_link_budget, que es el
% punto de paso universal del resto de pipelines).
assert_cfg_coherent(cfg);

verbose = getf(pt, 'verbose', false);
label   = getf(pt, 'label',   '');
keepCdf = getf(pt, 'keep_cdf', false);
elMinV  = getf(pt, 'elMin_valid', []);

%% 0b. Stream de aleatoriedad PROPIO DEL PUNTO (ver cabecera)
%  Threefry soporta substreams: el punto i usa el substream i de la misma semilla,
%  luego su resultado no depende de que worker lo ejecute ni en que orden.
seed  = getf(pt, 'seed', 42);
idxPt = getf(pt, 'idx',  1);
rs = RandStream('Threefry', 'Seed', seed);
rs.Substream = idxPt;
RandStream.setGlobalStream(rs);

%% 0c. Esquemas a evaluar en este punto
if isfield(pt,'cases') && ~isempty(pt.cases)
    cases = pt.cases;
    if ~iscell(cases), cases = {cases}; end
else
    one = struct();
    for f = {'scheme','Delta','alpha','adaptive','name'}
        if isfield(pt, f{1}), one.(f{1}) = pt.(f{1}); end
    end
    cases = {one};
end
nC = numel(cases);

% cfg de cada esquema, resuelta UNA vez (no depende del bloque temporal). Se parte
% SIEMPRE de la cfg del punto: los esquemas no se contaminan entre si aunque uno
% deje campos sin tocar.
cfgCase = cell(1,nC);
for c = 1:nC
    cc   = cases{c};
    cfgC = cfg;
    if isfield(cc,'scheme')   && ~isempty(cc.scheme),   cfgC.ffr.scheme   = cc.scheme;   end
    if isfield(cc,'Delta')    && ~isempty(cc.Delta),    cfgC.ffr.Delta    = cc.Delta;    end
    if isfield(cc,'alpha')    && ~isempty(cc.alpha) && ~any(isnan(cc.alpha))
        cfgC.ffr.alpha = cc.alpha;
    end
    if isfield(cc,'adaptive') && ~isempty(cc.adaptive), cfgC.ffr.adaptive = cc.adaptive; end
    if isfield(cc,'alpha_min')&& ~isempty(cc.alpha_min),cfgC.ffr.alpha_min= cc.alpha_min;end
    if isfield(cc,'alpha_max')&& ~isempty(cc.alpha_max),cfgC.ffr.alpha_max= cc.alpha_max;end
    cfgCase{c} = cfgC;
end

tPoint = tic;

%% 1. Constelacion, usuarios y layout (NO dependen del instante)
tvec  = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
Nt    = numel(tvec);
sats  = build_constellation(cfg);
users = build_user_grid(cfg);
BL    = build_beam_layout(cfg);
M     = users.M;  N = sats.N;  nB = BL.nBeams;

% Satelites ELEGIBLES como servidores. Con cfg.geom.serving_cid = [] (defecto) son
% todos y el comportamiento es el de siempre; en escenarios MULTI-OPERADOR (E4) se
% fija a 1 para que el usuario se sirva SOLO de su propio operador y los ajenos
% queden como interferentes puros (ver associate_serving).
servMask = [];
if isfield(cfg,'geom') && isfield(cfg.geom,'serving_cid') && ~isempty(cfg.geom.serving_cid)
    servMask = ismember(sats.cid, cfg.geom.serving_cid);
    if ~any(servMask)
        error('run_one_density:servingCid', ...
            'cfg.geom.serving_cid = [%s] no selecciona ningun satelite.', ...
            num2str(cfg.geom.serving_cid));
    end
end

% Tamano de bloque temporal (ver cabecera). nBlk = Nt -> un solo bloque = sin trocear.
[nBlk, nBlocks] = resolve_time_block(cfg, Nt);

%% 2. Acumuladores de MUESTRAS [M x Nt]
%  Se acumulan las MATRICES DE MUESTRA, no KPIs parciales. Ocupan M*Nt (unos
%  cientos de kB), despreciable frente al pico de M*N*Nt que se quiere acotar, y a
%  cambio compute_kpis recibe al final EXACTAMENTE el mismo array que recibiria sin
%  trocear: percentiles, CDF, Jain, R_agg y el veredicto de viabilidad salen
%  identicos sin necesidad de ninguna formula de combinacion por bloques (que es
%  justo donde es facil equivocarse: un p5 NO es la media de los p5 de cada bloque).
acc = repmat(struct('SINR',nan(M,Nt), 'CN',nan(M,Nt), 'pen',nan(M,Nt), ...
                    'Bu',nan(M,Nt), 'isC',false(M,Nt), ...
                    'alpha',nan(1,Nt), 'tau',nan(1,Nt), 'n0',nan(1,Nt)), 1, nC);
elevAcc  = nan(M,Nt);
thetaAcc = nan(M,Nt);
nVisAcc  = nan(M,Nt);
covAcc   = false(M,Nt);
metaC    = cell(1,nC);

%% 3. Bucle de BLOQUES TEMPORALES (pipeline dependiente del tiempo)
for b = 1:nBlocks
    kIdx = ((b-1)*nBlk + 1) : min(b*nBlk, Nt);
    tv   = tvec(kIdx);

    R_eci  = propagate(cfg, sats, tv);
    R_ecef = eci2ecef(cfg, R_eci, tv);
    % sats.minElev = mascara de elevacion POR CONSTELACION (Starlink 25 deg vs
    % OneWeb 55 deg). Con una sola constelacion equivale al escalar global y el
    % resultado es bit-identico al de antes.
    G      = compute_geometry(cfg, users, R_ecef, sats.minElev);
    S      = associate_serving(cfg, G, servMask);
    L      = compute_link_budget(cfg, S);
    INTF   = compute_interference(cfg, sats, users, G, S, L);
    ctx    = ffr_context(cfg, sats, users, BL, R_ecef, G, S, L, INTF, verbose && b == 1);

    % Liberar AQUI es lo que hace que el troceado ACOTE el pico y no solo lo
    % reparta: si G siguiera viva, el bloque siguiente sumaria su geometria a la
    % anterior y el maximo volveria a ser el de la ventana completa.
    clear G R_eci R_ecef INTF S L;

    for c = 1:nC
        A = ffr_allocate(cfgCase{c}, ctx, false);     % verbose = false
        F = compute_sinr_ffr(cfgCase{c}, ctx, A);

        acc(c).SINR(:,kIdx) = F.SINR_dB;
        acc(c).CN(:,kIdx)   = F.CN_dB;
        acc(c).pen(:,kIdx)  = F.penalty_dB;
        acc(c).Bu(:,kIdx)   = A.B_user_Hz;
        acc(c).isC(:,kIdx)  = A.isCenter;
        acc(c).alpha(kIdx)  = A.alpha;
        acc(c).tau(kIdx)    = A.tau;
        % CARGA MEDIDA de la celda (usuarios cubiertos / celdas ocupadas, por
        % instante). Es la n0 de la convencion sched='share' de ffr_allocate:
        % B_user = B_sub/(n0*f_clase). Se propaga porque es el denominador que
        % convierte ancho de SUB-BANDA en ancho POR USUARIO, y sin el no se puede
        % auditar un resultado expresado en tasa por usuario (lo necesita EB, que
        % deriva "cuantos usuarios caben"). Es GEOMETRIA pura: no depende del
        % esquema ni del ancho de banda. Campo ANADIDO, puramente informativo: no
        % entra en ningun calculo y no cambia ningun numero.
        acc(c).n0(kIdx)     = A.load_n0;

        % La SINR de referencia en reuso-1 (clasificador 'sinr' y tau por cuantil)
        % se calcula UNA vez POR BLOQUE y se cachea en ctx: es independiente del
        % esquema por construccion (ffr_allocate la obtiene con una asignacion de
        % reuso pleno) y se define por instante, luego cachearla dentro del bloque
        % da el mismo valor que sin trocear y garantiza que todos los esquemas
        % comparten EXACTAMENTE la misma clasificacion centro/borde.
        if isempty(ctx.SINR_ref_dB), ctx.SINR_ref_dB = A.SINR_ref_dB; end

        if b == 1
            metaC{c} = struct('Delta', A.Delta, 'sched', A.sched);
        end
    end

    elevAcc(:,kIdx)  = ctx.elev_deg;
    thetaAcc(:,kIdx) = ctx.theta_deg;
    nVisAcc(:,kIdx)  = ctx.nVis;
    covAcc(:,kIdx)   = ctx.cov;
    clear ctx;
end

%% 4. KPIs sobre las muestras COMPLETAS (una sola llamada, como sin trocear)
K  = cell(1,nC);
KV = cell(1,nC);
maskV = [];
if ~isempty(elMinV), maskV = covAcc & elevAcc >= elMinV; end

isCenter_ref = [];
casesOut = cell(1,nC);

for c = 1:nC
    FF = struct('SINR_dB', acc(c).SINR, 'CN_dB', acc(c).CN, 'penalty_dB', acc(c).pen);
    AL = struct('B_user_Hz', acc(c).Bu, 'isCenter', acc(c).isC, ...
                'scheme', cfgCase{c}.ffr.scheme, 'Delta', metaC{c}.Delta, ...
                'alpha', acc(c).alpha);

    Kc = compute_kpis(cfgCase{c}, FF, AL);
    if ~keepCdf, Kc = strip_cdf(Kc); end
    K{c} = Kc;

    if ~isempty(maskV)
        Kv = compute_kpis(cfgCase{c}, FF, AL, maskV);
        if ~keepCdf, Kv = strip_cdf(Kv); end
        KV{c} = Kv;
    end

    % Guarda: la poblacion de BORDE (sobre la que se juzga la viabilidad de H1)
    % debe ser la misma en todos los esquemas del punto.
    if isempty(isCenter_ref)
        isCenter_ref = acc(c).isC;
    elseif ~isequal(acc(c).isC, isCenter_ref)
        error('run_one_density:edgePop', ...
            ['La clasificacion centro/borde NO coincide entre esquemas en el punto ' ...
             '"%s": los KPIs de borde compararian poblaciones distintas.'], label);
    end

    % Definicion del caso, ya resuelta (para que OUT sea autoexplicativo)
    co = struct();
    co.name     = getf(cases{c}, 'name', ...
                       sprintf('%s(D=%d)', cfgCase{c}.ffr.scheme, cfgCase{c}.ffr.Delta));
    co.scheme   = cfgCase{c}.ffr.scheme;
    co.Delta    = metaC{c}.Delta;
    co.adaptive = cfgCase{c}.ffr.adaptive;
    co.alpha_cfg   = cfgCase{c}.ffr.alpha;
    co.alpha_mean  = mean(acc(c).alpha);
    co.alpha_range = [min(acc(c).alpha), max(acc(c).alpha)];
    co.tau_mean    = mean(acc(c).tau);
    co.sched       = metaC{c}.sched;
    co.load_n0_mean  = mean(acc(c).n0, 'omitnan');          % usuarios/celda (media temporal)
    co.load_n0_range = [min(acc(c).n0), max(acc(c).n0)];
    casesOut{c}    = co;
end

%% 5. Diagnostico del punto
cov = covAcc;
d = struct();
d.covFrac      = sum(cov(:)) / max(numel(cov),1);
if any(cov(:))
    ev = elevAcc(cov);
    d.elev_mean = mean(ev);  d.elev_min = min(ev);  d.elev_max = max(ev);
    % Offset angular del usuario a su centro de haz. La MEDIA mide la COMPRESION
    % ANGULAR del cluster Earth-fixed (a baja elevacion el cluster se ve mas
    % pequeno desde el satelite y los offsets encogen); el MAXIMO es la senal de
    % si la rejilla de usuarios sobresale del cluster (comparar con s = spacing_deg).
    d.theta_mean_deg = mean(thetaAcc(cov));
    d.theta_max_deg  = max(thetaAcc(cov));
    d.elev_p5   = prctile(ev, 5);
else
    d.elev_mean = NaN;  d.elev_min = NaN;  d.elev_max = NaN;
    d.theta_mean_deg = NaN;  d.theta_max_deg = NaN;  d.elev_p5 = NaN;
end
d.nVis_mean    = mean(nVisAcc(:), 'omitnan');
d.spacing_deg  = BL.spacing_deg;
d.spacing_km   = BL.spacing_km;
if ~isempty(maskV), d.nValid = sum(maskV(:)); else, d.nValid = NaN; end

% CRITERIO DE H2 (experimento E5): la sub-banda INTERIOR solo compensa su
% ancho si  Delta * SE_centro / SE_borde > 1.  Se registra en cada punto de
% densidad porque es justo lo que E3 debe vigilar: si al saturar cruza 1, H2
% recuperaria margen y el resultado negativo de E5 seria especifico de T=66.
%
% SOLO SE DEVUELVE PARA LOS ESQUEMAS 'ffr'; en el resto queda NaN. El criterio se
% deduce derivando  R = Bc*SE_c + Be*SE_b  respecto de alpha, luego presupone que
% existe una SUB-BANDA INTERIOR distinta de la de borde. En 'reuse1' y 'reuseD' no
% la hay: centro y borde ven el MISMO conjunto co-canal, el cociente sale > 1 por
% construccion y no significa nada. Devolverlo para todos los esquemas ya causo
% una lectura equivocada (el runner de OneWeb imprimio un "interior rentable: SI"
% espurio), asi que el campo se renombra a
% crit_ratio_ffr y se pone a NaN fuera de su dominio de validez: si alguien lo lee
% donde no aplica, obtiene NaN en vez de un numero plausible pero sin sentido.
d.crit_ratio_ffr = nan(1,nC);
for c = 1:nC
    if ~strcmpi(casesOut{c}.scheme, 'ffr'), continue; end     % fuera de dominio
    sec = K{c}.center.SE_mean_bpsHz;
    seb = K{c}.edge.SE_mean_bpsHz;
    if ~isnan(sec) && ~isnan(seb) && seb > 0
        d.crit_ratio_ffr(c) = casesOut{c}.Delta * sec / seb;
    end
end

%% 6. Salida
OUT = struct();
OUT.label   = label;
OUT.T       = cfg.constellations(1).T;
OUT.P       = cfg.constellations(1).P;
OUT.N       = sats.N;
OUT.M       = users.M;
OUT.Nt      = numel(tvec);
OUT.nBeams  = BL.nBeams;
OUT.elMin_valid = elMinV;
OUT.cases   = casesOut;
OUT.K       = K;
OUT.KV      = KV;
OUT.diag    = d;

% Troceado temporal y coste de memoria (informativo; NO afecta a los KPIs)
OUT.mem = struct( ...
    'timeBlock',     nBlk, ...
    'nBlocks',       nBlocks, ...
    'peak_model_GB', mem_peak_GB(M, N, nBlk, nB), ...   % con el troceado aplicado
    'peak_full_GB',  mem_peak_GB(M, N, Nt,   nB), ...   % si se hiciera de una pieza
    'samples_MB',    nC*5*M*Nt*8 / 2^20);               % coste de los acumuladores

% Metadatos NO deterministas: van APARTE para que la comparacion
% secuencial-vs-parfor pueda excluirlos sin excluir nada sustantivo.
OUT.runinfo = struct('elapsed_s', toc(tPoint), 'idx', idxPt, 'seed', seed);

if verbose
    fprintf('[punto %s] T=%d N=%d M=%d Nt=%d | cobertura %.2f | elev %.1f-%.1f deg | %.1f s\n', ...
        label, OUT.T, OUT.N, OUT.M, OUT.Nt, d.covFrac, d.elev_min, d.elev_max, ...
        OUT.runinfo.elapsed_s);
end
end

% -------------------------------------------------------------------------
function [nBlk, nBlocks] = resolve_time_block(cfg, Nt)
%RESOLVE_TIME_BLOCK  Instantes por bloque a partir de cfg.compute.timeBlock.
%   [] / Inf / <=0 / >= Nt  ->  un solo bloque (sin trocear, comportamiento por
%   defecto). Lectura defensiva: una cfg anterior a este campo (p.ej. recuperada
%   de un .mat antiguo) sigue funcionando sin trocear.
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
end

% -------------------------------------------------------------------------
function gb = mem_peak_GB(M, N, Ntb, nB)
%MEM_PEAK_GB  Pico de memoria de UN bloque temporal de Ntb instantes [GB].
%   Delega en mem_peak_model_GB, que es la FUENTE UNICA del modelo (74 B/elem
%   M*N*Ntb + 40 B/elem M*nBeams*Ntb, FASE A, cota superior). El mismo fichero lo
%   usa estimate_sweep_memory para dimensionar el barrido antes de lanzarlo, de
%   modo que lo estimado y lo medido no pueden divergir.
gb = mem_peak_model_GB(M, N, Ntb, nB);
end

% -------------------------------------------------------------------------
function K = strip_cdf(K)
%STRIP_CDF  Elimina las CDF de los KPIs (vectores largos) para aligerar la
%   transferencia worker -> cliente. No afecta a ninguna metrica: las CDF se
%   reconstruyen reejecutando el punto con pt.keep_cdf = true.
f = {'all','center','edge'};
for i = 1:numel(f)
    if isfield(K, f{i})
        if isfield(K.(f{i}), 'cdf_SINR'), K.(f{i}).cdf_SINR = []; end
        if isfield(K.(f{i}), 'cdf_R'),    K.(f{i}).cdf_R    = []; end
    end
end
end

% -------------------------------------------------------------------------
function v = getf(s, name, defaultVal)
%GETF  Campo opcional de un struct con valor por defecto.
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultVal;
end
end
