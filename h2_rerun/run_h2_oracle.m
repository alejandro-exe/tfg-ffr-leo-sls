%% RUN_H2_ORACLE  Cota ORACULO del reparto de banda: el mejor alpha(t) posible.
%
%  VERSION 2 (CORRECCION DE LA BUSQUEDA). Misma fisica, mismo modelo afin y misma
%  puerta anclada a los .mat que la v1; lo unico que cambia es COMO se busca.
%  Resultados en h2_oracle_v2.mat; h2_oracle.mat (v1) se CONSERVA como registro.
%
%  QUE ESTABA MAL EN LA v1, y es la causa raiz de los sep < 0:
%  el conjunto de candidatos de cada coordenada era [rejilla, a(t)+refinamientos] y
%  NO CONTENIA a(t). Es decir, cada coordenada estaba OBLIGADA a moverse aunque su
%  valor actual fuese el mejor, luego el ascenso NO era monotono y podia DEGRADAR su
%  propia semilla. Por eso en tau=60 el "oraculo" (10.175236) quedaba por debajo de
%  la adaptativa nominal (10.191244), que es factible y estaba entre las semillas.
%  Mas reinicios no lo habrian arreglado: no era falta de exploracion.
%
%  CORRECCIONES DE LA v2:
%   (a) a(t) entra en el conjunto de candidatos -> el ascenso es MONOTONO;
%   (b) SEMILLAS GARANTIZADAS en TODAS las optimizaciones (cabecera y los 9 niveles
%       de frontera): los alpha CONSTANTES de la rejilla refinada, la adaptativa
%       nominal y ademas el alpha constante que alcanza EXACTAMENTE el nivel
%       (bisección sobre el modelo afin). El resultado devuelto es el MAXIMO entre
%       lo encontrado y la mejor semilla FACTIBLE evaluada tal cual;
%   (c) 200 reinicios en el punto de cabecera y 40 por nivel de frontera, con
%       n_feasible, n_at_best e histograma de los factibles;
%   (d) dispersion de corr(alpha, elevMean) sobre TODOS los reinicios factibles a
%       menos del 2% del mejor objetivo, no solo la del mejor;
%   (e) FRONTERA ESTATICA REFINADA: barrido de alpha con paso 0.01 en
%       [alpha*-0.10, alpha*+0.10], con GR y matThr recalculados sobre ella. La
%       frontera es CONCAVA y la interpolacion por cuerdas la SUBESTIMA, luego con
%       la rejilla gruesa la ganancia sale inflada; ademas matThr es medio paso de
%       la rejilla y por tanto DEPENDE de ella. Se reportan las dos rejillas.
%
%  CONSECUENCIA EXIGIBLE de (a)+(b): front.sep >= 0 en todos los niveles. Si sale
%  negativo el script ABORTA: con la mejor semilla factible garantizada y la
%  frontera concava, un negativo solo puede ser un fallo de la busqueda.
%
%  MODELO AFIN (sin cambios respecto de la v1). Con sched='share', PSD constante y
%  tau fijo por cuantil: la SINR no depende de alpha y la clasificacion la fija tau,
%  luego  R_u(t) = k_u(t)*alpha(t) (centro) y k_u(t)*(1-alpha(t)) (borde), con
%  k_u(t) independiente de alpha. Basta UNA asignacion del motor para extraer
%  k_u(t); a partir de ahi cualquier alpha(t) se evalua en aritmetica pura.
%  Validado y BLOQUEANTE contra compute_kpis en los tres tau.
%
%  NO genera figuras. NO modifica nada fuera de h2_rerun/. NO toca el motor.
%
%  Uso: run_h2_oracle   (desde la raiz del proyecto o desde h2_rerun/).
%
%  TFG UC3M - Viabilidad de FFR adaptativa en constelaciones LEO.

clear; clc;

%% 0. Rutas y log -----------------------------------------------------------
H2DIR   = fileparts(mfilename('fullpath'));
ROOTDIR = fileparts(H2DIR);
RESDIR  = fullfile(H2DIR, 'results');
LOGDIR  = fullfile(H2DIR, 'logs');
if ~exist(RESDIR,'dir'), mkdir(RESDIR); end
if ~exist(LOGDIR,'dir'), mkdir(LOGDIR); end
addpath(ROOTDIR);

LOGF = fullfile(LOGDIR, 'run_h2_oracle_v2.txt');
diary off;
if exist(LOGF,'file'), delete(LOGF); end
diary(LOGF);
tRun = tic;

%% 1. Perfil: IDENTICO al de run_e5_adaptive / run_h2_anchors
cfg = config_default();
cfg.constellations(1).T = 1584;
cfg.constellations(1).P = 72;
cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 3;
cfg.time.dt       = 60;
cfg.time.duration = 7200;
cfg.ffr.tau_mode = 'quantile';
cfg.ffr.tau_q    = 50;
cfg.beams.nRings = 5;
cfg.compute.timeBlock = 10;

DELTA       = 3;
elMin_valid = 45;
alphas      = 0.05:0.05:0.95;       % rejilla GRUESA (la de E5; ancla de regresion)
TAUS        = [50 40 60];
ABOX        = [0.05 0.95];                                          % <-- DECISION
NLEV        = 9;                    % niveles de la frontera oraculo  <-- DECISION
NREST_HEAD  = 200;                  % reinicios en el punto de cabecera
NREST_FRONT = 40;                   % reinicios por nivel de frontera
REF_WIN     = 0.10;                 % semiventana del refinado de alpha  <-- DECISION
REF_STEP    = 0.01;                 % paso del refinado                  <-- DECISION
CORR_TOL    = 0.02;                 % 2% del mejor objetivo para la dispersion
SEED        = 20260831;                                             % <-- DECISION

fprintf('===== H2-ORACLE v2 | cota oraculo del reparto de banda =====\n');
fprintf('Ku DL: f=%.0f GHz | B=%.0f MHz | HPBW=%.2f deg | T=%d/P=%d | Delta=%d\n', ...
    cfg.radio.freq_GHz, cfg.radio.B_MHz, cfg.radio.beamwidth3dB_deg, ...
    cfg.constellations(1).T, cfg.constellations(1).P, DELTA);
fprintf('caja alpha = [%.2f, %.2f] | tau_q = %s | niveles = %d\n', ...
    ABOX(1), ABOX(2), mat2str(TAUS), NLEV);
fprintf('reinicios: %d (cabecera) / %d (por nivel) | refinado: paso %.2f en +-%.2f de alpha*\n', ...
    NREST_HEAD, NREST_FRONT, REF_STEP, REF_WIN);

tW = tic;  atm_loss_dB(cfg, 45);  tWarm = toc(tW);
fprintf('[warmup] tabla P.618 en %.1f s\n', tWarm);
tExp = tic;

%% 2. Contexto compartido
MEM = estimate_sweep_memory(cfg, {struct('T', cfg.constellations(1).T)}, 1, ...
                            struct('verbose', false));
[ctx, sats, users, BL, BINFO] = build_ctx_blocked(cfg, true);
fprintf('Rejilla: %d usuarios | Nt=%d | N=%d satelites | %d haces\n', ...
    users.M, ctx.Nt, sats.N, BL.nBeams);

cfg.ffr.scheme = 'reuse1';  cfg.ffr.Delta = 1;  cfg.ffr.adaptive = false;
A_r1 = ffr_allocate(cfg, ctx, false);
ctx.SINR_ref_dB = A_r1.SINR_ref_dB;
maskV = ctx.cov & ctx.elev_deg >= elMin_valid;
Nt    = ctx.Nt;

PAR = struct('alphas',alphas,'DELTA',DELTA,'ABOX',ABOX,'NLEV',NLEV, ...
             'NREST_HEAD',NREST_HEAD,'NREST_FRONT',NREST_FRONT, ...
             'REF_WIN',REF_WIN,'REF_STEP',REF_STEP,'CORR_TOL',CORR_TOL);

%% 3. Bloques de tau
nT = numel(TAUS);
OB = cell(1, nT);
for it = 1:nT
    fprintf('\n#######################################################################\n');
    fprintf('#  BLOQUE tau_q = %d%%\n', TAUS(it));
    fprintf('#######################################################################\n');
    OB{it} = oracle_block(cfg, ctx, maskV, PAR, TAUS(it), SEED + it);
end

%% 4. Salidas agregadas con los nombres pedidos
alpha_oracle     = nan(nT, Nt);
gain_oracle      = nan(1, nT);
gain_oracle_ref  = nan(1, nT);
corr_oracle_elev = nan(1, nT);
for it = 1:nT
    alpha_oracle(it,:)   = OB{it}.alpha_oracle;
    gain_oracle(it)      = OB{it}.gain_oracle;
    gain_oracle_ref(it)  = OB{it}.gain_oracle_ref;
    corr_oracle_elev(it) = OB{it}.corr_oracle_elev;
end
restarts     = OB{1}.restarts;
affine_check = OB{1}.affine_check;
front_oracle = OB{1}.front;
corr_spread  = OB{1}.corr_spread;
for it = 2:nT
    restarts(it)     = OB{it}.restarts;
    affine_check(it) = OB{it}.affine_check;
    front_oracle(it) = OB{it}.front;
    corr_spread(it)  = OB{it}.corr_spread;
end

%% 5. Resumen numerico
fprintf('\n=======================================================================\n');
fprintf('  RESUMEN\n');
fprintf('=======================================================================\n');
fprintf('%6s %8s %11s %12s | %11s %10s %9s | %11s %10s %9s\n', ...
    'tau_q','alpha*','objetivo','R5b oraculo', ...
    'front.gruesa','ganancia','material','front.refin.','ganancia','material');
fprintf('%s\n', repmat('-',1,116));
for it = 1:nT
    o = OB{it};
    fprintf('%5d%% %8.2f %11.6f %12.6f | %11.6f %+10.6f %9s | %11.6f %+10.6f %9s\n', ...
        o.tau_q, o.alpha_star, o.target, o.p5e_oracle, ...
        o.p5e_frontier,     o.gain_oracle,     ternary(o.gain_oracle     > o.matThr,    'SI','no'), ...
        o.p5e_frontier_ref, o.gain_oracle_ref, ternary(o.gain_oracle_ref > o.matThr_ref,'SI','no'));
end
fprintf('%s\n', repmat('-',1,116));
fprintf('umbral de materialidad: rejilla gruesa %s | rejilla refinada %s\n', ...
    mat2str(round(cellfun(@(o) o.matThr,     OB),6)), ...
    mat2str(round(cellfun(@(o) o.matThr_ref, OB),6)));

fprintf('\n%6s %10s %12s %12s %12s %12s %11s\n', ...
    'tau_q','reinicios','factibles','mejor','peor(fact.)','desv.tipica','n en mejor');
fprintf('%s\n', repmat('-',1,80));
for it = 1:nT
    r = OB{it}.restarts;
    fprintf('%5d%% %10d %12d %12.6f %12.6f %12.3e %11d\n', ...
        OB{it}.tau_q, r.n, r.n_feasible, r.best, r.worst, r.std, r.n_at_best);
end
fprintf('%s\n', repmat('-',1,80));

fprintf('\n--- corr(alpha, elevMean): del mejor y DISPERSION sobre los factibles a <%.0f%% ---\n', ...
    100*CORR_TOL);
fprintf('%6s %10s %10s %10s %10s %10s %10s\n', ...
    'tau_q','del mejor','n en banda','media','desv.tip.','minimo','maximo');
fprintf('%s\n', repmat('-',1,74));
for it = 1:nT
    c = OB{it}.corr_spread;
    fprintf('%5d%% %10.4f %10d %10.4f %10.4f %10.4f %10.4f\n', ...
        OB{it}.tau_q, OB{it}.corr_oracle_elev, c.n, c.mean, c.std, c.min, c.max);
end
fprintf('%s\n', repmat('-',1,74));
fprintf('corr(alpha_nominal, elevMean) = %s (la regla de E5, por construccion)\n', ...
    mat2str(round(cellfun(@(o) o.corr_nominal_elev, OB),4)));

fprintf('\n--- frontera oraculo frente a frontera estatica REFINADA ---\n');
fprintf('%6s %12s %12s %12s %12s %10s %10s\n', ...
    'tau_q','sep. max','en R5g','sep. media','sep. min','umbral','n material');
fprintf('%s\n', repmat('-',1,86));
for it = 1:nT
    f = OB{it}.front;
    fprintf('%5d%% %12.6f %12.6f %12.6f %12.6f %10.6f %10d\n', ...
        OB{it}.tau_q, f.sep_max, f.sep_max_at, f.sep_mean, f.sep_min, ...
        OB{it}.matThr_ref, f.n_material);
end
fprintf('%s\n', repmat('-',1,86));

%% 6. PUERTA DE REGRESION (anclada a los .mat, sin constantes transcritas)
fprintf('\n#######################################################################\n');
fprintf('#  PUERTA DE REGRESION  (referencias leidas de los .mat, tolerancia 1e-9)\n');
fprintf('#######################################################################\n');
TOLREG = 1e-9;
E5F = fullfile(ROOTDIR, 'e5_adaptive_results.mat');
ANF = fullfile(RESDIR,  'h2_anchors.mat');
for f = {E5F, ANF}
    if ~exist(f{1},'file')
        diary off;
        error('run_h2_oracle:sinRef', 'No se encuentra %s: sin el no hay puerta.', f{1});
    end
end
E5R = load(E5F, 'sw','alpha_star','GR','ratio_all','alphas','DELTA','KK');
ANR = load(ANF, 'tau_sweep','tau_list','alphas','DELTA');
fprintf('referencias: %s\n             %s\n', E5F, ANF);

if ~isequal(E5R.alphas(:).', alphas(:).') || E5R.DELTA ~= DELTA || ...
   ~isequal(ANR.alphas(:).', alphas(:).')
    diary off;
    error('run_h2_oracle:ejeRef', 'El eje de alpha o Delta no coinciden con los .mat.');
end

REG = struct('name',{{}}, 'ref',[], 'got',[], 'dif',[], 'ok',[]);
swf = {'R_p5','R_agg','R_p5_edge','R_p5_V','R_agg_V','R_p5_edge_V'};

i50 = find(TAUS == 50, 1);
for i = 1:numel(swf)
    REG = reg_add_max(REG, sprintf('[E5] sw.%s (afin vs motor)', swf{i}), ...
        max(abs(OB{i50}.sw.(swf{i}) - E5R.sw.(swf{i}))));
end
REG = reg_add(REG, '[E5] alpha_star',  E5R.alpha_star,  OB{i50}.alpha_star);
REG = reg_add(REG, '[E5] GR.gain_y',   E5R.GR.gain_y,   OB{i50}.gain_nominal);
REG = reg_add(REG, '[E5] KK{3}.all.R_p5_Mbps (objetivo)', ...
                                       E5R.KK{3}.all.R_p5_Mbps, OB{i50}.target);
REG = reg_add(REG, '[E5] KK{3}.edge.R_p5_Mbps', ...
                                       E5R.KK{3}.edge.R_p5_Mbps, OB{i50}.p5e_nominal);

for it = 1:nT
    ia = find(ANR.tau_list == TAUS(it), 1);
    if isempty(ia)
        diary off;
        error('run_h2_oracle:tauRef', 'h2_anchors.mat no trae el bloque tau_q = %d.', TAUS(it));
    end
    S = ANR.tau_sweep(ia);
    for i = 1:numel(swf)
        REG = reg_add_max(REG, sprintf('[ANC tau=%d] sw.%s', TAUS(it), swf{i}), ...
            max(abs(OB{it}.sw.(swf{i}) - S.sw.(swf{i}))));
    end
    REG = reg_add(REG, sprintf('[ANC tau=%d] alpha_star', TAUS(it)), S.alpha_star, OB{it}.alpha_star);
    REG = reg_add(REG, sprintf('[ANC tau=%d] matThr',     TAUS(it)), S.matThr,     OB{it}.matThr);
    REG = reg_add(REG, sprintf('[ANC tau=%d] GR.gain_y',  TAUS(it)), S.gain_y,     OB{it}.gain_nominal);
    REG = reg_add(REG, sprintf('[ANC tau=%d] R_p5_adapt', TAUS(it)), S.R_p5_glob_adapt, OB{it}.target);
    REG = reg_add(REG, sprintf('[ANC tau=%d] R_p5_edge_adapt', TAUS(it)), ...
                                                                  S.R_p5_edge_adapt, OB{it}.p5e_nominal);
    REG = reg_add(REG, sprintf('[ANC tau=%d] centerFrac', TAUS(it)), S.centerFrac, OB{it}.centerFrac);
    REG = reg_add(REG, sprintf('[ANC tau=%d] n_edge',     TAUS(it)), S.n_edge,     OB{it}.n_edge);
end

for it = 1:nT
    ac = OB{it}.affine_check;
    REG = reg_add_max(REG, sprintf('[AFIN tau=%d] max|dif| KPIs (3 casos)', TAUS(it)), ac.max_dif);
    REG = reg_add_max(REG, sprintf('[P5RAPIDO tau=%d] max|dif| vs prctile', TAUS(it)), ac.max_dif_p5fast);
    % La rejilla REFINADA debe reproducir la GRUESA en los alpha comunes: si no,
    % el refinado estaria cambiando el barrido en vez de densificarlo.
    REG = reg_add_max(REG, sprintf('[REFIN tau=%d] comun grueso-refinado', TAUS(it)), ...
        OB{it}.refine_consistency);
    % Exigencia de la correccion: la separacion no puede ser negativa en ningun nivel.
    REG = reg_add_max(REG, sprintf('[SEP>=0 tau=%d] min(front.sep) negativo', TAUS(it)), ...
        max(0, -min(OB{it}.front.sep)));
end

fprintf('\n%-46s %16s %16s %11s %s\n','magnitud','referencia','obtenido','|dif|','');
fprintf('%s\n', repmat('-',1,112));
nBad = 0;
for i = 1:numel(REG.name)
    d = abs(REG.got(i) - REG.ref(i));
    REG.dif(i) = d;  REG.ok(i) = d <= TOLREG;
    if ~REG.ok(i), nBad = nBad + 1; end
    fprintf('%-46s %16.9f %16.9f %11.2e %s\n', REG.name{i}, REG.ref(i), REG.got(i), d, ...
        ternary(REG.ok(i),'OK','*** FALLA ***'));
end
fprintf('%s\n', repmat('-',1,112));
REG.tol = TOLREG;  REG.nBad = nBad;
fprintf('%d de %d dentro de %.0e (max|dif| = %.2e).\n', ...
    sum(REG.ok), numel(REG.ok), TOLREG, max(REG.dif));

fprintf('\n-----------------------------------------------------------------------\n');
if nBad > 0
    fprintf('PUERTA FALLIDA: %d comprobacion(es) fuera de tolerancia. NO se guarda.\n', nBad);
    diary off;
    error('run_h2_oracle:regresion', ...
        '%d comprobacion(es) fuera de %.0e. Ver %s.', nBad, TOLREG, LOGF);
end
fprintf('PUERTA SUPERADA.\n');
fprintf('-----------------------------------------------------------------------\n');

%% 7. Guardado
timing = struct('warmup_p618_s', tWarm, 'experiment_s', toc(tExp), 'total_s', toc(tRun));
fprintf('\n[tiempos] experimento %.1f s (+ %.1f s de warmup) | total %.1f s\n', ...
    timing.experiment_s, timing.warmup_p618_s, timing.total_s);

OUTF = fullfile(RESDIR, 'h2_oracle_v2.mat');
save(OUTF, 'alpha_oracle','gain_oracle','gain_oracle_ref','restarts','corr_oracle_elev', ...
     'corr_spread','affine_check','front_oracle','OB','TAUS','alphas','DELTA','ABOX', ...
     'PAR','SEED','cfg','BL','REG','timing','MEM','BINFO');
fprintf('\nResultados en %s\n', OUTF);
fprintf('(h2_oracle.mat, la tanda v1, se CONSERVA sin tocar)\n');
fprintf('Log en %s\n', LOGF);
fprintf('\nFIN.\n');
diary off;

%% ========================= funciones locales =========================

function O = oracle_block(cfg, ctx, maskV, PAR, tau_q, seed)
%ORACLE_BLOCK  Validacion afin + frontera refinada + cota oraculo para UN tau.

cfg.ffr.tau_mode = 'quantile';
cfg.ffr.tau_q    = tau_q;
cfg.ffr.scheme   = 'ffr';
cfg.ffr.Delta    = PAR.DELTA;
Nt     = ctx.Nt;
alphas = PAR.alphas;
DELTA  = PAR.DELTA;
ABOX   = PAR.ABOX;

%% (1) Llamadas REALES al motor: extraccion + los tres casos de validacion
cfg.ffr.adaptive = false;
cfg.ffr.alpha    = 0.50;                 % alpha DEDICADO a la extraccion de k
A0 = ffr_allocate(cfg, ctx, false);
F0 = compute_sinr_ffr(cfg, ctx, A0);

cfg.ffr.alpha = 0.55;
A55 = ffr_allocate(cfg, ctx, false);  F55 = compute_sinr_ffr(cfg, ctx, A55);
K55 = compute_kpis(cfg, F55, A55);
cfg.ffr.alpha = 0.40;
A40 = ffr_allocate(cfg, ctx, false);  F40 = compute_sinr_ffr(cfg, ctx, A40);
K40 = compute_kpis(cfg, F40, A40);
cfg.ffr.adaptive = true;  cfg.ffr.adapt_dir = +1;
AAD = ffr_allocate(cfg, ctx, false);  FAD = compute_sinr_ffr(cfg, ctx, AAD);
KAD = compute_kpis(cfg, FAD, AAD);
cfg.ffr.adaptive = false;

isC   = A0.isCenter;
cov   = ctx.cov;
elevM = AAD.POL.elevMean_deg;
alpha_nominal = AAD.alpha;

%% (2) Extraccion de k_u(t)
R0   = A0.B_user_Hz .* log2(1 + 10.^(F0.SINR_dB/10));
a0   = A0.alpha;
Kmat = nan(size(R0));
A0m  = repmat(a0, size(R0,1), 1);
Kmat( isC & cov) = R0( isC & cov) ./ A0m( isC & cov);
Kmat(~isC & cov) = R0(~isC & cov) ./ (1 - A0m(~isC & cov));

P.G  = pack_mask(cov,          Kmat, isC, Nt);
P.E  = pack_mask(cov & ~isC,   Kmat, isC, Nt);
P.GV = pack_mask(maskV,        Kmat, isC, Nt);
P.EV = pack_mask(maskV & ~isC, Kmat, isC, Nt);
P.aggG  = agg_prep(cov,   Kmat, isC, Nt);
P.aggGV = agg_prep(maskV, Kmat, isC, Nt);

%% (3) VALIDACION DEL MODELO AFIN
cases = {'estatica alpha=0.55', 0.55*ones(1,Nt), K55; ...
         'estatica alpha=0.40', 0.40*ones(1,Nt), K40; ...
         'adaptativa nominal',  alpha_nominal,   KAD};
fprintf('\n--- validacion del modelo afin (k extraido con alpha=0.50) ---\n');
fprintf('%-22s %14s %14s %14s\n','caso','|d R_p5 glob|','|d R_p5 borde|','|d R_agg|');
fprintf('%s\n', repmat('-',1,68));
dmax = 0;  dtab = nan(3,3);
for c = 1:3
    av = cases{c,2};  Kr = cases{c,3};
    [g, e, ag] = aff_kpis(P, av);
    dtab(c,:) = [abs(g - Kr.all.R_p5_Mbps), abs(e - Kr.edge.R_p5_Mbps), ...
                 abs(ag - Kr.all.R_agg_Mbps)];
    dmax = max(dmax, max(dtab(c,:)));
    fprintf('%-22s %14.2e %14.2e %14.2e\n', cases{c,1}, dtab(c,1), dtab(c,2), dtab(c,3));
end
fprintf('%s\nmax|dif| = %.3e\n', repmat('-',1,68), dmax);

%% (4) Verificacion del percentil rapido contra prctile (en Mbps)
dp5 = 0;
rs  = RandStream('threefry','Seed',seed);
for c = 1:6
    if c <= 3, av = cases{c,2}; else, av = ABOX(1) + diff(ABOX)*rand(rs,1,Nt); end
    vG = aff_values(P.G, av) / 1e6;   vE = aff_values(P.E, av) / 1e6;
    dp5 = max(dp5, abs(p5_sorted(sort(vG), numel(vG), 5) - prctile(vG,5)));
    dp5 = max(dp5, abs(p5_sorted(sort(vE), numel(vE), 5) - prctile(vE,5)));
end
fprintf('percentil rapido vs prctile: max|dif| = %.3e (6 secuencias alpha)\n', dp5);
O.affine_check = struct('max_dif', dmax, 'per_case', dtab, ...
    'cases', {{cases{1,1}, cases{2,1}, cases{3,1}}}, 'max_dif_p5fast', dp5, ...
    'alpha_extraccion', 0.50);

%% (5) Barrido estatico GRUESO (ancla de regresion) y REFINADO
sw  = affine_sweep(P, alphas);
[~, iBest] = max(sw.R_p5);
alpha_star = alphas(iBest);

% Rejilla REFINADA: la gruesa + paso REF_STEP en [alpha*-REF_WIN, alpha*+REF_WIN].
aR = (alpha_star - PAR.REF_WIN) : PAR.REF_STEP : (alpha_star + PAR.REF_WIN);
alphasR = unique(round([alphas, aR]*1e6)/1e6);
alphasR = alphasR(alphasR >= ABOX(1) - 1e-12 & alphasR <= ABOX(2) + 1e-12);
swR = affine_sweep(P, alphasR);
[~, iBestR] = max(swR.R_p5);
alpha_starR = alphasR(iBestR);

% Consistencia: en los alpha COMUNES la refinada debe reproducir la gruesa.
[~, ia, ib] = intersect(round(alphas*1e6), round(alphasR*1e6));
refine_consistency = max(abs(sw.R_p5_edge(ia) - swR.R_p5_edge(ib)));
refine_consistency = max(refine_consistency, max(abs(sw.R_p5(ia) - swR.R_p5(ib))));

fprintf('\nrejilla gruesa: %d alpha | rejilla refinada: %d alpha (paso %.2f en [%.2f, %.2f])\n', ...
    numel(alphas), numel(alphasR), PAR.REF_STEP, min(aR), max(aR));
fprintf('alpha* grueso = %.2f | alpha* refinado = %.2f | consistencia en comunes = %.2e\n', ...
    alpha_star, alpha_starR, refine_consistency);

%% (6) Fronteras y umbrales de materialidad, en las DOS rejillas
[gN, eN] = aff_kpis(P, alpha_nominal);
GR      = frontier_gain(sw.R_p5,  sw.R_p5_edge,  gN, eN);
GRr     = frontier_gain(swR.R_p5, swR.R_p5_edge, gN, eN);
resolY  = median(abs(diff(sw.R_p5_edge)));
resolYr = median(abs(diff(swR.R_p5_edge)));
matThr     = max(0.02*abs(GR.y_frontier),  0.5*resolY);
matThr_ref = max(0.02*abs(GRr.y_frontier), 0.5*resolYr);

fprintf('\nadaptativa nominal: R_p5 = %.6f | R_p5_borde = %.6f\n', gN, eN);
fprintf('  frontera GRUESA   en ese R_p5 = %.6f | ganancia = %+.6f | umbral = %.6f\n', ...
    GR.y_frontier,  GR.gain_y,  matThr);
fprintf('  frontera REFINADA en ese R_p5 = %.6f | ganancia = %+.6f | umbral = %.6f\n', ...
    GRr.y_frontier, GRr.gain_y, matThr_ref);
fprintf('  la cuerda gruesa SUBESTIMA la frontera en %.6f Mbps (frontera concava)\n', ...
    GRr.y_frontier - GR.y_frontier);

%% (7) COTA ORACULO en el objetivo de la adaptativa nominal
target = gN;
OPT    = build_opt(P, ABOX);
fprintf('\n--- cota oraculo | objetivo R_p5 global >= %.6f | %d reinicios ---\n', ...
    target, PAR.NREST_HEAD);

starts = build_starts(P, OPT, alphasR, alpha_nominal, target, ABOX, ...
                      PAR.NREST_HEAD, alpha_starR, seed+1000);
[aOra, fOra, rinfo] = oracle_search(OPT, starts, target);

[gO, eO, agO] = aff_kpis(P, aOra);
GRo  = frontier_gain(sw.R_p5,  sw.R_p5_edge,  gO, eO);
GRor = frontier_gain(swR.R_p5, swR.R_p5_edge, gO, eO);
fprintf('reinicios = %d | factibles = %d | mejor = %.6f | peor(fact.) = %.6f | desv = %.3e | n en mejor = %d\n', ...
    rinfo.n, rinfo.n_feasible, rinfo.best, rinfo.worst, rinfo.std, rinfo.n_at_best);
print_hist(rinfo.f_feasible, 'objetivo de los reinicios FACTIBLES [Mbps]');
fprintf('alpha oraculo: R_p5 global = %.6f (objetivo %.6f) | R_p5 borde = %.6f | R_agg = %.4f\n', ...
    gO, target, eO, agO);
fprintf('  vs frontera GRUESA   %.6f -> ganancia %+.6f (umbral %.6f) %s\n', ...
    GRo.y_frontier,  GRo.gain_y,  matThr,     ternary(GRo.gain_y  > matThr,    'MATERIAL','no material'));
fprintf('  vs frontera REFINADA %.6f -> ganancia %+.6f (umbral %.6f) %s\n', ...
    GRor.y_frontier, GRor.gain_y, matThr_ref, ternary(GRor.gain_y > matThr_ref,'MATERIAL','no material'));
fprintf('rango de alpha oraculo = [%.4f, %.4f] | media = %.4f\n', ...
    min(aOra), max(aOra), mean(aOra));

%% (8) Correlacion con la elevacion: la del mejor y su DISPERSION
co = corr_nan(aOra, elevM);
cn = corr_nan(alpha_nominal, elevM);
band = rinfo.feasible & (rinfo.f_all >= rinfo.best*(1 - PAR.CORR_TOL));
cvals = nan(1, sum(band));
ib = find(band);
for q = 1:numel(ib), cvals(q) = corr_nan(rinfo.alpha_all(ib(q),:), elevM); end
cvals = cvals(isfinite(cvals));
corr_spread = struct('n', numel(cvals), 'mean', mean_or_nan(cvals), ...
    'std', std_or_nan(cvals), 'min', min_or_nan(cvals), 'max', max_or_nan(cvals), ...
    'values', cvals, 'tol', PAR.CORR_TOL);
fprintf(['corr(alpha_oraculo, elevMean) = %+.4f | corr(alpha_nominal, elevMean) = %+.4f\n' ...
         'dispersion sobre los %d reinicios factibles a <%.0f%% del mejor: ' ...
         'media %+.4f | desv %.4f | rango [%+.4f, %+.4f]\n'], ...
    co, cn, corr_spread.n, 100*PAR.CORR_TOL, corr_spread.mean, corr_spread.std, ...
    corr_spread.min, corr_spread.max);

%% (9) FRONTERA ORACULO sobre la rejilla REFINADA
lev = linspace(min(swR.R_p5), max(swR.R_p5), PAR.NLEV);
NL  = PAR.NLEV;
front = struct('levels', lev, 'p5e', nan(1,NL), 'p5g', nan(1,NL), ...
               'stat', nan(1,NL), 'sep', nan(1,NL), 'stat_coarse', nan(1,NL), ...
               'sep_coarse', nan(1,NL), 'nrest', PAR.NREST_FRONT, ...
               'alpha', nan(NL,Nt), 'n_feasible', nan(1,NL), 'n_at_best', nan(1,NL), ...
               'rbest', nan(1,NL), 'rworst', nan(1,NL));
fprintf('\n--- frontera oraculo (%d niveles, %d reinicios por nivel, rejilla refinada) ---\n', ...
    NL, PAR.NREST_FRONT);
fprintf('%12s %12s %12s %12s %12s %8s %8s %8s\n', ...
    'objetivo','R5g logrado','R5b oraculo','R5b estatica','separacion','fact.','n mejor','material');
fprintf('%s\n', repmat('-',1,96));
for L = 1:NL
    st = build_starts(P, OPT, alphasR, alpha_nominal, lev(L), ABOX, ...
                      PAR.NREST_FRONT, alpha_starR, seed + 2000 + L);
    st{end+1} = aOra; %#ok<AGROW>   arranque templado desde el punto de cabecera
    [aL, ~, rL] = oracle_search(OPT, st, lev(L));
    [gL, eL] = aff_kpis(P, aL);
    front.p5e(L) = eL;  front.p5g(L) = gL;
    front.stat(L)        = interp_clamped(GRr.fx, GRr.fy, gL);
    front.sep(L)         = eL - front.stat(L);
    front.stat_coarse(L) = interp_clamped(GR.fx, GR.fy, gL);
    front.sep_coarse(L)  = eL - front.stat_coarse(L);
    front.alpha(L,:)     = aL;
    front.n_feasible(L)  = rL.n_feasible;
    front.n_at_best(L)   = rL.n_at_best;
    front.rbest(L)  = rL.best;  front.rworst(L) = rL.worst;
    fprintf('%12.6f %12.6f %12.6f %12.6f %+12.6f %8d %8d %8s\n', ...
        lev(L), gL, eL, front.stat(L), front.sep(L), rL.n_feasible, rL.n_at_best, ...
        ternary(front.sep(L) > matThr_ref,'SI','no'));
end
fprintf('%s\n', repmat('-',1,96));
[front.sep_max, iS] = max(front.sep);
front.sep_max_at = front.p5g(iS);
front.sep_mean   = mean(front.sep);
front.sep_min    = min(front.sep);
front.n_material = sum(front.sep > matThr_ref);
fprintf('separacion maxima = %+.6f en R5g = %.6f | media = %+.6f | minima = %+.6f\n', ...
    front.sep_max, front.sep_max_at, front.sep_mean, front.sep_min);
fprintf('niveles con separacion > umbral (%.6f): %d de %d\n', matThr_ref, front.n_material, NL);
fprintf('separacion media contra la frontera GRUESA (para contraste): %+.6f\n', ...
    mean(front.sep_coarse));

%% (10) Salida del bloque
O.tau_q         = tau_q;
O.sw            = sw;            O.swR = swR;
O.alphasR       = alphasR;
O.alpha_star    = alpha_star;    O.alpha_starR = alpha_starR;
O.refine_consistency = refine_consistency;
O.matThr        = matThr;        O.matThr_ref = matThr_ref;
O.target        = target;
O.p5e_nominal   = eN;
O.gain_nominal  = GR.gain_y;     O.gain_nominal_ref = GRr.gain_y;
O.alpha_oracle  = aOra;
O.p5g_oracle    = gO;            O.p5e_oracle = eO;
O.p5e_frontier  = GRo.y_frontier;    O.gain_oracle     = GRo.gain_y;
O.p5e_frontier_ref = GRor.y_frontier; O.gain_oracle_ref = GRor.gain_y;
O.f_oracle      = fOra;
O.restarts      = rinfo;
O.corr_oracle_elev  = co;
O.corr_nominal_elev = cn;
O.corr_spread   = corr_spread;
O.alpha_nominal = alpha_nominal;
O.elevMean      = elevM;
O.front         = front;
O.centerFrac    = K55.centerFrac;
O.n_edge        = K55.edge.n;
O.GRstat        = struct('fx', GR.fx,  'fy', GR.fy);
O.GRstatR       = struct('fx', GRr.fx, 'fy', GRr.fy);
end

% -------------------------------------------------------------------------
function sw = affine_sweep(P, av)
%AFFINE_SWEEP  Barrido de alpha CONSTANTE en aritmetica afin (los 6 campos de sw).
n = numel(av);
sw = struct('alpha', av, 'R_p5', nan(1,n), 'R_agg', nan(1,n), 'R_p5_edge', nan(1,n), ...
            'R_p5_V', nan(1,n), 'R_agg_V', nan(1,n), 'R_p5_edge_V', nan(1,n));
Nt = P.G.Nt;
for i = 1:n
    a = av(i)*ones(1,Nt);
    [g, e, ag]    = aff_kpis(P, a);
    [gv, ev, agv] = aff_kpis_V(P, a);
    sw.R_p5(i) = g;  sw.R_p5_edge(i) = e;  sw.R_agg(i) = ag;
    sw.R_p5_V(i) = gv;  sw.R_p5_edge_V(i) = ev;  sw.R_agg_V(i) = agv;
end
end

% -------------------------------------------------------------------------
function st = build_starts(P, OPT, alphasR, alpha_nominal, target, ABOX, nWant, aStarR, seed)
%BUILD_STARTS  SEMILLAS GARANTIZADAS de una optimizacion.
%   Siempre entran: (1) TODOS los alpha constantes de la rejilla refinada -- que son
%   los vertices con los que se construye la frontera estatica, luego el oraculo no
%   puede quedar por debajo de ninguno de ellos; (2) la adaptativa nominal; (3) el
%   alpha CONSTANTE que alcanza EXACTAMENTE el nivel objetivo, obtenido por biseccion
%   sobre el modelo afin. La (3) es la que hace exigible sep >= 0: da un punto de la
%   curva estatica VERDADERA en x = target, y como la curva es concava la cuerda que
%   interpola los vertices queda por debajo de ella.
%   El resto, hasta nWant, son aleatorios reproducibles.
Nt = P.G.Nt;
st = {};
for i = 1:numel(alphasR), st{end+1} = alphasR(i)*ones(1,Nt); end %#ok<AGROW>
st{end+1} = alpha_nominal;
aL = alpha_for_level(P, target, ABOX(1), aStarR);
st{end+1} = aL*ones(1,Nt);
rs = RandStream('threefry','Seed',seed);
while numel(st) < nWant
    st{end+1} = ABOX(1) + diff(ABOX)*rand(rs,1,Nt); %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function a = alpha_for_level(P, L, lo, hi)
%ALPHA_FOR_LEVEL  alpha CONSTANTE con R_p5 global = L, por biseccion.
%   R_p5 global crece con alpha en [0.05, alpha*], que es el tramo de la frontera.
%   Se devuelve el extremo FACTIBLE del intervalo final (p5g >= L).
if p5g_const(P, lo) >= L, a = lo; return; end
if p5g_const(P, hi) <= L, a = hi; return; end
for it = 1:40
    m = 0.5*(lo + hi);
    if p5g_const(P, m) < L, lo = m; else, hi = m; end
end
a = hi;
end

% -------------------------------------------------------------------------
function g = p5g_const(P, a)
g = prctile(aff_values(P.G, a*ones(1,P.G.Nt)) / 1e6, 5);
end

% -------------------------------------------------------------------------
function S = pack_mask(mask, Kmat, isC, Nt)
%PACK_MASK  Aplana las muestras de una mascara agrupadas POR INSTANTE (contiguas).
idx   = find(mask);
S.k   = Kmat(idx);
S.isC = isC(idx);
S.n   = numel(idx);
cnt   = sum(mask, 1);
S.ie  = cumsum(cnt(:)).';
S.is  = S.ie - cnt(:).' + 1;
S.cnt = cnt;  S.Nt = Nt;
S.t   = zeros(S.n,1);
for k = 1:Nt
    if cnt(k) > 0, S.t(S.is(k):S.ie(k)) = k; end
end
end

% -------------------------------------------------------------------------
function A = agg_prep(mask, Kmat, isC, Nt)
A.Sc = zeros(1,Nt);  A.Se = zeros(1,Nt);
for k = 1:Nt
    mk = mask(:,k);
    if ~any(mk), continue; end
    A.Sc(k) = sum(Kmat(mk &  isC(:,k), k));
    A.Se(k) = sum(Kmat(mk & ~isC(:,k), k));
end
A.active = sum(mask,1) > 0;
end

% -------------------------------------------------------------------------
function v = aff_values(S, alpha)
a = alpha(S.t).';
v = S.k .* (S.isC .* a + (~S.isC) .* (1 - a));
end

% -------------------------------------------------------------------------
function [p5g, p5e, ragg] = aff_kpis(P, alpha)
p5g  = prctile(aff_values(P.G, alpha) / 1e6, 5);
p5e  = prctile(aff_values(P.E, alpha) / 1e6, 5);
s    = alpha .* P.aggG.Sc + (1 - alpha) .* P.aggG.Se;
ragg = mean(s(P.aggG.active)) / 1e6;
end

% -------------------------------------------------------------------------
function [p5g, p5e, ragg] = aff_kpis_V(P, alpha)
p5g  = prctile(aff_values(P.GV, alpha) / 1e6, 5);
p5e  = prctile(aff_values(P.EV, alpha) / 1e6, 5);
s    = alpha .* P.aggGV.Sc + (1 - alpha) .* P.aggGV.Se;
ragg = mean(s(P.aggGV.active)) / 1e6;
end

% -------------------------------------------------------------------------
function v = p5_sorted(xs, n, p)
%P5_SORTED  Percentil p de n valores dado un PREFIJO ORDENADO xs (los menores).
%   Reproduce la definicion de prctile y se VERIFICA contra ella sobre los datos
%   reales antes de usarse.
h = n*p/100 + 0.5;
if h <= 1, v = xs(1); return; end
i = floor(h);  f = h - i;
if i + 1 > numel(xs)
    if numel(xs) >= n, v = xs(end); return; end
    error('p5_sorted:prefijo', ...
        ['El prefijo ordenado tiene %d elementos y el percentil %g de %d necesita ' ...
         'el indice %d: la poda es demasiado agresiva.'], numel(xs), p, n, i+1);
end
v = xs(i) + f*(xs(i+1) - xs(i));
end

% -------------------------------------------------------------------------
function OPT = build_opt(P, ABOX)
OPT.G = P.G;  OPT.E = P.E;
OPT.box = ABOX;
OPT.KG = ceil(P.G.n*0.05 + 0.5) + 2;
OPT.KE = ceil(P.E.n*0.05 + 0.5) + 2;
OPT.Nt = P.G.Nt;
OPT.grid = linspace(ABOX(1), ABOX(2), 10);                          % <-- DECISION
OPT.ref  = [-0.06 -0.03 -0.015 0.015 0.03 0.06];                    % <-- DECISION
OPT.maxPass = 5;                                                    % <-- DECISION
OPT.tol     = 1e-7;
OPT.thFac   = 1.5;      % holgura del umbral de poda (con respaldo EXACTO)
end

% -------------------------------------------------------------------------
function [aBest, fBest, info] = oracle_search(OPT, starts, target)
%ORACLE_SEARCH  Ascenso por coordenadas con reinicios y SEMILLAS GARANTIZADAS.
%   Devuelve el MAXIMO entre lo encontrado por el ascenso y la mejor SEMILLA
%   FACTIBLE evaluada tal cual: con el ascenso ya monotono esto es redundante, y
%   por eso mismo se deja -- es la red que impide devolver algo peor que un punto
%   que se tenia de partida (que es lo que fallaba en la v1).
nR   = numel(starts);
Nt   = OPT.Nt;
fAll = -inf(1,nR);  gAll = -inf(1,nR);  aAll = nan(nR,Nt);
fSeed = -inf(1,nR); gSeed = -inf(1,nR);
for r = 1:nR
    s0 = min(max(starts{r}(:).', OPT.box(1)), OPT.box(2));
    [fSeed(r), gSeed(r)] = obj_of(aff_values(OPT.G, s0), aff_values(OPT.E, s0), ...
                                  OPT.G.n, OPT.E.n);
    [a, f, g] = coord_ascent(OPT, s0, target);
    % Red de seguridad: quedarse con la semilla si el ascenso no la mejoro.
    if better(fSeed(r), gSeed(r), f, g, target)
        a = s0;  f = fSeed(r);  g = gSeed(r);
    end
    aAll(r,:) = a;  fAll(r) = f;  gAll(r) = g;
end
feas = gAll >= target - 1e-12;
if any(feas)
    fSel = fAll;  fSel(~feas) = -inf;
    [fBest, ib] = max(fSel);
else
    [~, ib] = max(gAll);  fBest = fAll(ib);
end
aBest = aAll(ib,:);
ff = fAll(feas);
info = struct('n', nR, 'n_feasible', sum(feas), 'best', fBest, ...
    'worst', min_or_nan(ff), 'std', std_or_nan(ff), 'mean', mean_or_nan(ff), ...
    'n_at_best', sum(abs(ff - fBest) <= 1e-9), ...
    'f_all', fAll, 'g_all', gAll, 'feasible', feas, 'alpha_all', aAll, ...
    'f_feasible', ff, 'f_seed', fSeed, 'g_seed', gSeed, ...
    'n_seed_wins', sum(fSeed >= fAll - 1e-12 & gSeed >= target - 1e-12));
end

% -------------------------------------------------------------------------
function [a, fVal, gVal] = coord_ascent(OPT, a0, target)
%COORD_ASCENT  Ascenso por coordenadas MONOTONO hasta convergencia.
%   CORRECCION v2: a(t) forma parte del conjunto de candidatos, luego una
%   coordenada solo se mueve si MEJORA. En la v1 no estaba y el ascenso podia
%   empeorar su propia semilla.
a  = min(max(a0(:).', OPT.box(1)), OPT.box(2));
G  = OPT.G;  E = OPT.E;
vG = aff_values(G, a);  vE = aff_values(E, a);
[fVal, gVal] = obj_of(vG, vE, G.n, E.n);

% Umbrales de poda iniciales (se adaptan solos; el respaldo mantiene la EXACTITUD)
sg = sort(vG);  thG = sg(min(OPT.KG, numel(sg))) * OPT.thFac;
se = sort(vE);  thE = se(min(OPT.KE, numel(se))) * OPT.thFac;

for pass = 1:OPT.maxPass
    f0 = fVal;  g0 = gVal;
    for t = 1:OPT.Nt
        gs = G.is(t);  ge = G.ie(t);
        es = E.is(t);  ee = E.ie(t);
        if ge < gs, continue; end

        % --- parte FIJA podada a los KG/KE menores. La poda es EXACTA: si el
        % umbral no deja bastantes elementos, se cae al orden completo.
        [Ag, thG] = pruned_prefix(vG, gs, ge, OPT.KG, thG, OPT.thFac);
        if ee >= es
            [Ae, thE] = pruned_prefix(vE, es, ee, OPT.KE, thE, OPT.thFac);
        else
            [Ae, thE] = pruned_prefix(vE, 1, 0, OPT.KE, thE, OPT.thFac);
        end

        kG = G.k(gs:ge);  cG = G.isC(gs:ge);
        kE = E.k(es:ee);

        cand = [OPT.grid, a(t), a(t) + OPT.ref];        % <-- a(t) INCLUIDO (v2)
        cand = unique(min(max(cand, OPT.box(1)), OPT.box(2)));
        bf = -inf;  bg = -inf;  ba = a(t);
        for c = 1:numel(cand)
            aa = cand(c);
            xg = sort([Ag; kG .* (cG*aa + (~cG)*(1-aa))]);
            gg = p5_sorted(xg, G.n, 5) / 1e6;
            if ee >= es
                xe = sort([Ae; kE * (1-aa)]);
            else
                xe = Ae;
            end
            ff = p5_sorted(xe, E.n, 5) / 1e6;
            if better(ff, gg, bf, bg, target), bf = ff; bg = gg; ba = aa; end
        end
        if ba ~= a(t)
            a(t) = ba;
            vG(gs:ge) = kG .* (cG*ba + (~cG)*(1-ba));
            if ee >= es, vE(es:ee) = kE * (1-ba); end
        end
        fVal = bf;  gVal = bg;
    end
    if (fVal - f0) <= OPT.tol && (gVal - g0) <= OPT.tol, break; end
end
end

% -------------------------------------------------------------------------
function [A, th] = pruned_prefix(v, is, ie, K, th, fac)
%PRUNED_PREFIX  Los K menores de v EXCLUYENDO el bloque [is,ie].
%   Acelera con un umbral adaptativo pero NO aproxima: si el umbral no deja al
%   menos K elementos, se ordena el vector completo. El percentil que se calcula
%   despues es por tanto exacto.
sel = v <= th;
if ie >= is, sel(is:ie) = false; end
A = sort(v(sel));
if numel(A) < K
    w = v;
    if ie >= is, w(is:ie) = []; end
    A = sort(w);
end
A = A(1:min(K, numel(A)));
if numel(A) >= K, th = A(K) * fac; end
end

% -------------------------------------------------------------------------
function tf = better(f, g, fRef, gRef, target)
%BETTER  Comparacion LEXICOGRAFICA factible-primero.
okN = g >= target;  okR = gRef >= target;
if okN && ~okR,  tf = true;  return; end
if ~okN && okR,  tf = false; return; end
if okN,          tf = f > fRef;  return; end
tf = g > gRef;
end

% -------------------------------------------------------------------------
function [f, g] = obj_of(vG, vE, nG, nE)
g = p5_sorted(sort(vG), nG, 5) / 1e6;
f = p5_sorted(sort(vE), nE, 5) / 1e6;
end

% -------------------------------------------------------------------------
function print_hist(x, titulo)
%PRINT_HIST  Histograma de texto (10 clases) de los objetivos factibles.
x = x(isfinite(x));
fprintf('  histograma | %s  (n = %d)\n', titulo, numel(x));
if isempty(x), fprintf('    (vacio)\n'); return; end
if range_of(x) <= 0
    fprintf('    %10.6f | %s %d\n', x(1), repmat('#', 1, min(50,numel(x))), numel(x));
    return;
end
nb = 10;
ed = linspace(min(x), max(x), nb+1);
cn = histcounts(x, ed);
mx = max(cn);
for b = 1:nb
    fprintf('    [%9.6f, %9.6f) %s %d\n', ed(b), ed(b+1), ...
        repmat('#', 1, round(40*cn(b)/max(mx,1))), cn(b));
end
end

% -------------------------------------------------------------------------
function FR = frontier_gain(x, y, xq, yq)
%FRONTIER_GAIN  COPIA de la funcion local de run_e5_adaptive.m (no invocable desde
%   fuera). Anclada por la puerta: gain_nominal sale de aqui.
x = x(:).';  y = y(:).';
ok = isfinite(x) & isfinite(y);
x = x(ok);   y = y(ok);
keep = true(1, numel(x));
for i = 1:numel(x)
    dom = (x >= x(i)) & (y >= y(i)) & ((x > x(i)) | (y > y(i)));
    if any(dom), keep(i) = false; end
end
fx = x(keep);  fy = y(keep);
[fx, is] = sort(fx);  fy = fy(is);
[fxu, ~, gidx] = unique(fx);
fyu = accumarray(gidx(:), fy(:), [], @max).';
FR.fx = fxu;  FR.fy = fyu;
FR.x_q = xq;  FR.y_q = yq;
FR.y_frontier = interp_clamped(fxu, fyu, xq);
FR.x_frontier = interp_clamped(fliplr(fyu), fliplr(fxu), yq);
FR.gain_y = yq - FR.y_frontier;
FR.gain_x = xq - FR.x_frontier;
FR.gain_y_rel = FR.gain_y / max(abs(FR.y_frontier), eps);
FR.gain_x_rel = FR.gain_x / max(abs(FR.x_frontier), eps);
FR.dominated = any((x >= xq) & (y >= yq) & ((x > xq) | (y > yq)));
FR.above     = ~FR.dominated && (FR.gain_y > 0);
end

% -------------------------------------------------------------------------
function yq = interp_clamped(xx, yy, xq)
if numel(xx) < 2
    if isempty(xx), yq = NaN; else, yq = yy(1); end
    return;
end
if xq <= xx(1),   yq = yy(1);   return; end
if xq >= xx(end), yq = yy(end); return; end
yq = interp1(xx, yy, xq, 'linear');
end

% -------------------------------------------------------------------------
function r = corr_nan(a, b)
ok = isfinite(a) & isfinite(b);
if sum(ok) < 3 || std(a(ok)) == 0 || std(b(ok)) == 0, r = NaN; return; end
c = corrcoef(a(ok), b(ok));  r = c(1,2);
end

% -------------------------------------------------------------------------
function REG = reg_add(REG, name, ref, got)
REG.name{end+1} = name;  REG.ref(end+1) = ref;  REG.got(end+1) = got;
REG.dif(end+1) = NaN;    REG.ok(end+1)  = false;
end

% -------------------------------------------------------------------------
function REG = reg_add_max(REG, name, maxdif)
REG.name{end+1} = name;  REG.ref(end+1) = 0;  REG.got(end+1) = maxdif;
REG.dif(end+1) = NaN;    REG.ok(end+1)  = false;
end

% -------------------------------------------------------------------------
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function v = min_or_nan(x),  if isempty(x), v = NaN; else, v = min(x);  end, end
function v = max_or_nan(x),  if isempty(x), v = NaN; else, v = max(x);  end, end
function v = mean_or_nan(x), if isempty(x), v = NaN; else, v = mean(x); end, end
function v = std_or_nan(x),  if isempty(x), v = NaN; else, v = std(x);  end, end
function v = range_of(x),    v = max(x) - min(x); end
