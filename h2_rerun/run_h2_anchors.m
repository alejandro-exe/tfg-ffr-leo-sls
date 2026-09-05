%% RUN_H2_ANCHORS  E5 reejecutado con las ANCLAS de la regla adaptativa liberadas.
%
%  QUE PREGUNTA RESPONDE. La regla adaptativa de E5 interpola alpha entre dos
%  ANCLAS DE ELEVACION fijadas en cfg a 25 y 90 grados:
%
%      w(t)     = clip((elev(t) - elMin_deg)/(elMax_deg - elMin_deg), 0, 1)
%      alpha(t) = alpha_min + (alpha_max - alpha_min)*w(t)          [adapt_dir=+1]
%
%  Pero con la constelacion REAL (T=1584) la elevacion MEDIA del servidor solo
%  recorre 52.23-85.53 deg, luego w(t) se queda en [0.419, 0.931] y alpha(t) en
%  [0.4676, 0.6725] de un nominal [0.3, 0.7]: la regla usa el 51% de su rango.
%  Ademas, el barrido de familia de anclas de E5 mueve (alpha_min, alpha_max) con
%  las anclas de ELEVACION CONGELADAS en 25/90, asi que nunca se ha comprobado si
%  el resultado negativo de H2 depende de tener la regla ESTRANGULADA.
%
%  Este runner lo comprueba. NO es un experimento nuevo: es E5 con la misma
%  configuracion exacta, sin figuras, mas cuatro bloques de diagnostico:
%    (A) dos casos con las anclas de elevacion AJUSTADAS al recorrido real
%        (elMin = min(elevMean), elMax = max(elevMean)), en los dos sentidos;
%    (B) barrido 2D de anclas de elevacion, 4x4 = 16 combinaciones;
%    (C) barrido de anclas de ALPHA con las anclas de elevacion ya ajustadas;
%    (D) barrido ESTATICO de tau (tau_q = 40/50/60), repitiendo el bloque completo.
%
%  QUE **NO** HACE:
%    - no genera figuras (ni figure, ni print, ni exportgraphics);
%    - no toca el motor: ffr_policy, ffr_allocate, compute_sinr_ffr y compute_kpis
%      se usan tal cual y todo se varia por cfg;
%    - no adapta tau: tau es ESTATICO en los tres bloques del punto (D). Lo unico
%      que cambia entre ellos es el CUANTIL fijo del reparto centro/borde;
%    - no usa K.alpha, que devuelve alpha(1) (limitacion B14 documentada) y es
%      ENGANOSO en modo adaptativo. Los estadisticos de alpha salen de ALLOC.alpha.
%
%  PUERTA DE REGRESION (obligatoria, antes de guardar): las cuatro configuraciones
%  originales y el barrido de alpha deben reproducir e5_adaptive_results.mat con
%  |dif| <= 1e-9. Si alguna falla, aborta y NO guarda. Es lo que garantiza que los
%  bloques nuevos se miden sobre el MISMO escenario publicado.
%
%  NOTA DE DEUDA TECNICA (declarada, como en run_convergence_nt_step4): las
%  funciones locales frontier_gain / interp_clamped / corr_nan / ternary son COPIA
%  de las de run_e5_adaptive.m, que al ser locales de un script no son invocables
%  desde fuera. La copia queda ANCLADA por la puerta de regresion: GR.gain_y,
%  GRi.gain_y y scan.bestP5 dependen de frontier_gain, luego si la copia divergiera
%  el script aborta.
%
%  Uso: run_h2_anchors   (desde la raiz del proyecto o desde h2_rerun/).
%  Resultados en h2_rerun/results/h2_anchors.mat | log en h2_rerun/logs/.
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
addpath(ROOTDIR);                      % el motor vive en la raiz del proyecto

LOGF = fullfile(LOGDIR, 'run_h2_anchors.txt');
diary off;
if exist(LOGF,'file'), delete(LOGF); end
diary(LOGF);

tRun = tic;

%% 1. Perfil de estudio: IDENTICO al de run_e5_adaptive (comparabilidad directa)
%  Copiado literalmente. Cualquier desviacion la cazaria la puerta de regresion.
cfg = config_default();

cfg.constellations(1).T = 1584;        % Walker 53:1584/72/1 (Starlink Shell-1)
cfg.constellations(1).P = 72;

cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 3;

cfg.time.dt       = 60;                % s
cfg.time.duration = 7200;              % s (2 h)

cfg.ffr.tau_mode = 'quantile';         % reparto centro/borde fijo por cuantil
cfg.ffr.tau_q    = 50;                 % nominal (el punto (D) barre 40/50/60)

cfg.beams.nRings = 5;                  % 91 celdas Earth-fixed
cfg.compute.timeBlock = 10;            % obligatorio a densidad real

DELTA       = 3;
elMin_valid = 45;                      % deg  submuestra de "geometria valida"
alphas      = 0.05:0.05:0.95;          % barrido fino: traza la FRONTERA estatica

% Anclas NOMINALES, capturadas ANTES de que nada las toque. Son el denominador de
% la "fraccion del rango nominal usada" y la referencia de los barridos.
ALPHA_MIN_NOM = cfg.ffr.alpha_min;     % 0.3
ALPHA_MAX_NOM = cfg.ffr.alpha_max;     % 0.7
EL_MIN_NOM    = cfg.ffr.elMin_deg;     % 25
EL_MAX_NOM    = cfg.ffr.elMax_deg;     % 90

% Coeficiente de la recta R_p5_borde(alpha) = a*(1-alpha), publicado en E5 y usado
% para el "alpha fijo EQUIVALENTE" de una configuracion adaptativa:
%     alpha_equiv = 1 - R_p5_borde / A_EDGE
% Se verifica mas abajo contra el ajuste sobre los 19 puntos del barrido.
% SOLO VALE PARA tau_q = 50: con otro cuantil cambia la poblacion de borde y por
% tanto la pendiente, asi que en el punto (D) se reporta la pendiente de cada bloque.
A_EDGE = 16.511205;                                                % <-- publicado

fprintf('===== H2-RERUN | ANCLAS DE LA REGLA ADAPTATIVA (E5 sin figuras) =====\n');
fprintf('Ku DL: f=%.0f GHz | B=%.0f MHz | EIRPdens=%.1f dBW/MHz | HPBW=%.2f deg | G/T=%.1f dB/K\n', ...
    cfg.radio.freq_GHz, cfg.radio.B_MHz, cfg.radio.EIRPdensity_dBWMHz, ...
    cfg.radio.beamwidth3dB_deg, cfg.radio.GT_dBK);
fprintf('Anclas NOMINALES: alpha [%.2f, %.2f] | elevacion [%.0f, %.0f] deg | adapt_dir=%+d\n', ...
    ALPHA_MIN_NOM, ALPHA_MAX_NOM, EL_MIN_NOM, EL_MAX_NOM, cfg.ffr.adapt_dir);

tW = tic;  atm_loss_dB(cfg, 45);  tWarm = toc(tW);
fprintf('[warmup] tabla P.618 construida en %.1f s (coste fijo, se descuenta)\n', tWarm);
tExp = tic;

%% 2. Contexto compartido (una sola vez para TODO el script)
MEM = estimate_sweep_memory(cfg, {struct('T', cfg.constellations(1).T)}, 1, ...
                            struct('verbose', true));

[ctx, sats, users, BL, BINFO] = build_ctx_blocked(cfg, true);
tvec = ctx.tvec;

fprintf('Rejilla local: %d usuarios (radio %.0f km, paso %.0f km) | Nt=%d | N=%d satelites\n', ...
    users.M, cfg.ground.radius_km, cfg.ground.step_km, ctx.Nt, sats.N);
fprintf('Layout: %d haces | s = %.4f deg = %.2f km\n', BL.nBeams, BL.spacing_deg, BL.spacing_km);

% SINR de referencia (reuso-1): UNA vez, cacheada en ctx, para que TODOS los
% esquemas y TODOS los bloques de tau clasifiquen con la MISMA referencia.
cfg.ffr.scheme = 'reuse1';  cfg.ffr.Delta = 1;  cfg.ffr.adaptive = false;
A_r1 = ffr_allocate(cfg, ctx);
ctx.SINR_ref_dB = A_r1.SINR_ref_dB;

maskV  = ctx.cov & ctx.elev_deg >= elMin_valid;
nCov_k = sum(ctx.cov, 1);              % [1 x Nt] usuarios cubiertos por instante
nCovTot = sum(ctx.cov(:));

%% ========================================================================
%  3. BLOQUE NOMINAL (tau_q = 50): barrido de alpha + 4 configuraciones
%  ========================================================================
%  Se ejecuta a traves de la MISMA funcion local que usara el barrido de tau del
%  punto (D). Fuente unica: si el bloque nominal reproduce E5 (puerta de
%  regresion), los bloques de tau=40/60 salen del mismo camino de codigo.
fprintf('\n#######################################################################\n');
fprintf('#  BLOQUE NOMINAL  tau_q = %d%%\n', cfg.ffr.tau_q);
fprintf('#######################################################################\n');

[TB50, cfg] = tau_block(cfg, ctx, maskV, alphas, DELTA, 50, true);

sw         = TB50.sw;
alpha_star = TB50.alpha_star;
cases      = TB50.cases;
KK         = TB50.KK;    KV = TB50.KV;
AA         = TB50.AA;    FF = TB50.FF;
FR = TB50.FR;   FRv = TB50.FRv;   FRi = TB50.FRi;
GR = TB50.GR;   GRv = TB50.GRv;   GRi = TB50.GRi;
matThr     = TB50.matThr;
ratio_bin  = TB50.ratio_bin;
ratio_all  = TB50.ratio_all;
iR1 = 1;  iST = 2;  iAD = 3;  iIN = 4;

%% 4. VALIDACION: caso DEGENERADO alpha_min == alpha_max ------------------
%  Identica a E5. Si la adaptativa no degenera EXACTAMENTE en la estatica, hay un
%  bug que la favorece y nada de lo que sigue es interpretable.
fprintf('\n=============== VALIDACION: CASO DEGENERADO (alpha_min = alpha_max) ===============\n');
cfg_d = cfg;
cfg_d.ffr.scheme    = 'ffr';   cfg_d.ffr.Delta = DELTA;
cfg_d.ffr.adaptive  = true;
cfg_d.ffr.alpha_min = alpha_star;
cfg_d.ffr.alpha_max = alpha_star;
A_d = ffr_allocate(cfg_d, ctx, false);
F_d = compute_sinr_ffr(cfg_d, ctx, A_d);
K_d = compute_kpis(cfg_d, F_d, A_d);

dAlpha = max(abs(A_d.alpha - AA{iST}.alpha));
dSINR  = max(abs(F_d.SINR_dB(ctx.cov) - FF{iST}.SINR_dB(ctx.cov)));
dBu    = max(abs(A_d.B_user_Hz(ctx.cov) - AA{iST}.B_user_Hz(ctx.cov)));
dCls   = sum(A_d.isCenter(:) ~= AA{iST}.isCenter(:));
dKpi   = max(abs([K_d.all.R_p5_Mbps  - KK{iST}.all.R_p5_Mbps, ...
                  K_d.edge.R_p5_Mbps - KK{iST}.edge.R_p5_Mbps, ...
                  K_d.all.R_agg_Mbps - KK{iST}.all.R_agg_Mbps]));
degenOK = (dAlpha==0) && (dSINR==0) && (dBu==0) && (dCls==0) && (dKpi==0);
fprintf('  max|d alpha| = %.3g | max|d SINR| = %.3g dB | max|d B_user| = %.3g Hz\n', ...
    dAlpha, dSINR, dBu);
fprintf('  d clasificacion = %d muestras | max|d KPI| = %.3g Mbps  ==> %s\n', ...
    dCls, dKpi, ternary(degenOK,'PASA','FALLA'));

%% 5. BARRIDO ORIGINAL DE ANCLAS DE ALPHA (elevacion CONGELADA en 25/90) --
%  Es el barrido de E5, reproducido tal cual: entra en la puerta de regresion
%  (scan.bestP5) y es la LINEA BASE contra la que se compara el barrido del
%  punto (C), que repite lo mismo con las anclas de elevacion ya ajustadas.
fprintf('\n===== (base) FAMILIA DE ANCLAS DE ALPHA, elevacion congelada en %g/%g deg =====\n', ...
    EL_MIN_NOM, EL_MAX_NOM);
aLo = 0.05:0.10:0.55;
aHi = 0.05:0.10:0.95;
cfg_s = cfg;  cfg_s.ffr.scheme = 'ffr';  cfg_s.ffr.Delta = DELTA;
cfg_s.ffr.adaptive = true;  cfg_s.ffr.adapt_dir = +1;
scan = alpha_anchor_scan(cfg_s, ctx, sw, aLo, aHi);
fprintf('  MEJOR ganancia (R_agg, R_p5_borde): %+.6f Mbps  (anclas %.2f-%.2f)\n', ...
    scan.bestAgg, scan.bestAgg_anchors(1), scan.bestAgg_anchors(2));
fprintf('  MEJOR ganancia (R_p5,  R_p5_borde): %+.6f Mbps  (anclas %.2f-%.2f)\n', ...
    scan.bestP5, scan.bestP5_anchors(1), scan.bestP5_anchors(2));

%% ========================================================================
%  6. (A) ANCLAS DE ELEVACION AJUSTADAS AL RECORRIDO REAL
%  ========================================================================
POLa   = AA{iAD}.POL;
elev_t = POLa.elevMean_deg;
alp_t  = AA{iAD}.alpha;

elFitMin = min(elev_t);
elFitMax = max(elev_t);

fprintf('\n#######################################################################\n');
fprintf('#  (A) ANCLAS DE ELEVACION AJUSTADAS AL RECORRIDO REAL\n');
fprintf('#######################################################################\n');
fprintf(['DIAGNOSTICO DE LA REGLA ESTRANGULADA (ejecucion nominal, anclas %g/%g):\n' ...
         '  elevacion media del servidor : [%.4f, %.4f] deg\n' ...
         '  peso w(t) = clip(...)        : [%.4f, %.4f]  (deberia ser [0,1] si las\n' ...
         '                                  anclas cubrieran el recorrido real)\n' ...
         '  alpha(t)                     : [%.4f, %.4f] de un nominal [%.2f, %.2f]\n' ...
         '  fraccion del rango usada     : %.2f%%\n'], ...
    EL_MIN_NOM, EL_MAX_NOM, elFitMin, elFitMax, min(POLa.w), max(POLa.w), ...
    min(alp_t), max(alp_t), ALPHA_MIN_NOM, ALPHA_MAX_NOM, ...
    100*(max(alp_t)-min(alp_t))/(ALPHA_MAX_NOM-ALPHA_MIN_NOM));
fprintf('ANCLAS AJUSTADAS: elMin_deg = %.6f | elMax_deg = %.6f (min y max de elevMean)\n', ...
    elFitMin, elFitMax);

anchors_used = struct( ...
    'elMin_nom', EL_MIN_NOM, 'elMax_nom', EL_MAX_NOM, ...
    'alpha_min_nom', ALPHA_MIN_NOM, 'alpha_max_nom', ALPHA_MAX_NOM, ...
    'elev_min_real', elFitMin, 'elev_max_real', elFitMax, ...
    'elMin_fit', elFitMin, 'elMax_fit', elFitMax, ...
    'w_range_nom', [min(POLa.w) max(POLa.w)], ...
    'alpha_range_nom', [min(alp_t) max(alp_t)], ...
    'frac_range_nom', (max(alp_t)-min(alp_t))/(ALPHA_MAX_NOM-ALPHA_MIN_NOM));

% Los dos casos nuevos, evaluados contra la MISMA frontera estatica (sw) del
% bloque nominal: es lo que hace la comparacion legitima.
cfg_a = cfg;
cfg_a.ffr.scheme = 'ffr';  cfg_a.ffr.Delta = DELTA;  cfg_a.ffr.adaptive = true;
cfg_a.ffr.alpha_min = ALPHA_MIN_NOM;  cfg_a.ffr.alpha_max = ALPHA_MAX_NOM;
cfg_a.ffr.elMin_deg = elFitMin;       cfg_a.ffr.elMax_deg = elFitMax;

cfg_a.ffr.adapt_dir = +1;
A_p = ffr_allocate(cfg_a, ctx, false);
F_p = compute_sinr_ffr(cfg_a, ctx, A_p);
K_p = compute_kpis(cfg_a, F_p, A_p);
GR_anch = frontier_gain(sw.R_p5, sw.R_p5_edge, K_p.all.R_p5_Mbps, K_p.edge.R_p5_Mbps);
FR_anch = frontier_gain(sw.R_agg/1e3, sw.R_p5_edge, K_p.all.R_agg_Mbps/1e3, K_p.edge.R_p5_Mbps);

cfg_a.ffr.adapt_dir = -1;
A_m = ffr_allocate(cfg_a, ctx, false);
F_m = compute_sinr_ffr(cfg_a, ctx, A_m);
K_m = compute_kpis(cfg_a, F_m, A_m);
GRi_anch = frontier_gain(sw.R_p5, sw.R_p5_edge, K_m.all.R_p5_Mbps, K_m.edge.R_p5_Mbps);
FRi_anch = frontier_gain(sw.R_agg/1e3, sw.R_p5_edge, K_m.all.R_agg_Mbps/1e3, K_m.edge.R_p5_Mbps);

%% 7. DIAGNOSTICOS DE alpha(t) (los que E5 no reporta) --------------------
%  IMPORTANTE: NO se usa K.alpha, que devuelve alpha(1) (limitacion B14) y en modo
%  adaptativo es enganoso. Todo sale de ALLOC.alpha, que es el vector [1 x Nt] real.
cfgNames = {'adapt NOMINAL (25/90, +1)', 'adapt NOMINAL invertida (25/90, -1)', ...
            'adapt ANCLAS AJUSTADAS (+1)', 'adapt ANCLAS AJUSTADAS (-1)'};
alphaVecs = {AA{iAD}.alpha, AA{iIN}.alpha, A_p.alpha, A_m.alpha};
Kalist    = {KK{iAD}, KK{iIN}, K_p, K_m};
GAlist    = {GR, GRi, GR_anch, GRi_anch};

alpha_stats = struct('name',{{}}, 'mean_t',[], 'mean_samples',[], 'min',[], 'max',[], ...
                     'range',[], 'frac_nominal_range',[], 'alpha_equiv',[], ...
                     'R_p5_glob',[], 'R_p5_edge',[], 'R_agg_Gbps',[], 'gain_y',[]);
for c = 1:numel(alphaVecs)
    av = alphaVecs{c};
    alpha_stats.name{end+1}            = cfgNames{c};
    alpha_stats.mean_t(end+1)          = mean(av);
    alpha_stats.mean_samples(end+1)    = sum(av .* nCov_k) / nCovTot;
    alpha_stats.min(end+1)             = min(av);
    alpha_stats.max(end+1)             = max(av);
    alpha_stats.range(end+1)           = max(av) - min(av);
    alpha_stats.frac_nominal_range(end+1) = (max(av)-min(av)) / (ALPHA_MAX_NOM-ALPHA_MIN_NOM);
    alpha_stats.alpha_equiv(end+1)     = 1 - Kalist{c}.edge.R_p5_Mbps / A_EDGE;
    alpha_stats.R_p5_glob(end+1)       = Kalist{c}.all.R_p5_Mbps;
    alpha_stats.R_p5_edge(end+1)       = Kalist{c}.edge.R_p5_Mbps;
    alpha_stats.R_agg_Gbps(end+1)      = Kalist{c}.all.R_agg_Mbps/1e3;
    alpha_stats.gain_y(end+1)          = GAlist{c}.gain_y;
end
alpha_stats.A_EDGE = A_EDGE;
alpha_stats.matThr = matThr;

fprintf('\n--- DIAGNOSTICOS DE alpha(t)  (de ALLOC.alpha; K.alpha NO se usa, ver B14) ---\n');
fprintf('%-30s %9s %9s %19s %8s %8s %9s\n', ...
    'configuracion','media_t','media_mu','rango [min,max]','usado','a_equiv','R_5bord');
fprintf('%s\n', repmat('-',1,100));
for c = 1:numel(alphaVecs)
    fprintf('%-30s %9.4f %9.4f  [%.4f, %.4f] %7.1f%% %8.4f %9.3f\n', ...
        alpha_stats.name{c}, alpha_stats.mean_t(c), alpha_stats.mean_samples(c), ...
        alpha_stats.min(c), alpha_stats.max(c), 100*alpha_stats.frac_nominal_range(c), ...
        alpha_stats.alpha_equiv(c), alpha_stats.R_p5_edge(c));
end
fprintf('%s\n', repmat('-',1,100));
fprintf(['  media_t  = media de alpha(t) sobre INSTANTES | media_mu = media sobre MUESTRAS\n' ...
         '  cubiertas (ponderada por usuarios cubiertos por instante).\n' ...
         '  usado    = (max-min)/(alpha_max_nom - alpha_min_nom), fraccion del rango nominal.\n' ...
         '  a_equiv  = 1 - R_p5_borde/%.6f, el alpha FIJO que daria ese R_p5 de borde\n' ...
         '             (la recta R_borde(alpha) se verifica mas abajo).\n'], A_EDGE);

fprintf('\n--- (A) POSICION FRENTE A LA MISMA FRONTERA ESTATICA (umbral material %.6f) ---\n', matThr);
fprintf('%-30s %10s %11s %11s %11s %10s %s\n', ...
    'configuracion','R_p5_glob','R_p5_borde','frontera','ganancia','R_agg[Gb]','material?');
fprintf('%s\n', repmat('-',1,105));
for c = 1:numel(alphaVecs)
    G = GAlist{c};
    fprintf('%-30s %10.4f %11.4f %11.4f %+11.6f %10.4f  %s\n', ...
        alpha_stats.name{c}, alpha_stats.R_p5_glob(c), alpha_stats.R_p5_edge(c), ...
        G.y_frontier, G.gain_y, alpha_stats.R_agg_Gbps(c), ...
        ternary(G.above && G.gain_y > matThr, 'SI', 'no'));
end
fprintf('%s\n', repmat('-',1,105));
fprintf(['LECTURA: con las anclas AJUSTADAS la regla usa el 100%% de su rango de alpha\n' ...
         '  (recorre [%.2f, %.2f] en vez de [%.4f, %.4f]). Si aun asi no supera la frontera,\n' ...
         '  el resultado negativo de H2 NO se debia a tener la regla estrangulada.\n'], ...
    ALPHA_MIN_NOM, ALPHA_MAX_NOM, anchors_used.alpha_range_nom(1), anchors_used.alpha_range_nom(2));

%% ========================================================================
%  8. (B) BARRIDO 2D DE ANCLAS DE ELEVACION  (4 x 4 = 16)
%  ========================================================================
fprintf('\n#######################################################################\n');
fprintf('#  (B) BARRIDO 2D DE ANCLAS DE ELEVACION (alpha_min/alpha_max NOMINALES)\n');
fprintf('#######################################################################\n');
elLo = [25 40 50 55];
elHi = [70 80 86 90];
scan_el = struct('elLo',elLo, 'elHi',elHi, ...
    'gain_agg', nan(numel(elLo),numel(elHi)), 'gain_p5', nan(numel(elLo),numel(elHi)), ...
    'alpha_min',nan(numel(elLo),numel(elHi)), 'alpha_max',nan(numel(elLo),numel(elHi)), ...
    'frac_range',nan(numel(elLo),numel(elHi)), ...
    'R_p5',nan(numel(elLo),numel(elHi)), 'R_p5_edge',nan(numel(elLo),numel(elHi)), ...
    'R_agg',nan(numel(elLo),numel(elHi)), 'material',false(numel(elLo),numel(elHi)));

cfg_e = cfg;
cfg_e.ffr.scheme = 'ffr';  cfg_e.ffr.Delta = DELTA;  cfg_e.ffr.adaptive = true;
cfg_e.ffr.adapt_dir = +1;
cfg_e.ffr.alpha_min = ALPHA_MIN_NOM;  cfg_e.ffr.alpha_max = ALPHA_MAX_NOM;

fprintf('%6s %6s | %9s %9s %19s %8s | %11s %s\n', ...
    'elMin','elMax','R_p5','R_5borde','rango de alpha','usado','ganancia','material?');
fprintf('%s\n', repmat('-',1,96));
for i = 1:numel(elLo)
    for j = 1:numel(elHi)
        cfg_e.ffr.elMin_deg = elLo(i);
        cfg_e.ffr.elMax_deg = elHi(j);
        Ae = ffr_allocate(cfg_e, ctx, false);
        Fe = compute_sinr_ffr(cfg_e, ctx, Ae);
        Ke = compute_kpis(cfg_e, Fe, Ae);
        g1 = frontier_gain(sw.R_agg/1e3, sw.R_p5_edge, Ke.all.R_agg_Mbps/1e3, Ke.edge.R_p5_Mbps);
        g2 = frontier_gain(sw.R_p5,      sw.R_p5_edge, Ke.all.R_p5_Mbps,      Ke.edge.R_p5_Mbps);
        scan_el.gain_agg(i,j)  = g1.gain_y;
        scan_el.gain_p5(i,j)   = g2.gain_y;
        scan_el.alpha_min(i,j) = min(Ae.alpha);
        scan_el.alpha_max(i,j) = max(Ae.alpha);
        scan_el.frac_range(i,j)= (max(Ae.alpha)-min(Ae.alpha))/(ALPHA_MAX_NOM-ALPHA_MIN_NOM);
        scan_el.R_p5(i,j)      = Ke.all.R_p5_Mbps;
        scan_el.R_p5_edge(i,j) = Ke.edge.R_p5_Mbps;
        scan_el.R_agg(i,j)     = Ke.all.R_agg_Mbps;
        scan_el.material(i,j)  = g2.above && (g2.gain_y > matThr);
        fprintf('%6.0f %6.0f | %9.4f %9.4f  [%.4f, %.4f] %7.1f%% | %+11.6f  %s\n', ...
            elLo(i), elHi(j), Ke.all.R_p5_Mbps, Ke.edge.R_p5_Mbps, ...
            scan_el.alpha_min(i,j), scan_el.alpha_max(i,j), 100*scan_el.frac_range(i,j), ...
            g2.gain_y, ternary(scan_el.material(i,j),'SI','no'));
    end
end
fprintf('%s\n', repmat('-',1,96));
[scan_el.bestP5,  kE] = max(scan_el.gain_p5(:));
[iE,jE] = ind2sub(size(scan_el.gain_p5), kE);
scan_el.bestP5_anchors = [elLo(iE) elHi(jE)];
[scan_el.bestAgg, kA2] = max(scan_el.gain_agg(:));
[iA2,jA2] = ind2sub(size(scan_el.gain_agg), kA2);
scan_el.bestAgg_anchors = [elLo(iA2) elHi(jA2)];
scan_el.anyMaterial = any(scan_el.material(:));
fprintf('  MEJOR ganancia (R_p5, R_p5_borde) de las 16: %+.6f Mbps en (elMin=%g, elMax=%g)\n', ...
    scan_el.bestP5, scan_el.bestP5_anchors(1), scan_el.bestP5_anchors(2));
fprintf('  MEJOR ganancia (R_agg, R_p5_borde) de las 16: %+.6f Mbps en (elMin=%g, elMax=%g)\n', ...
    scan_el.bestAgg, scan_el.bestAgg_anchors(1), scan_el.bestAgg_anchors(2));
fprintf('  Alguna combinacion MATERIAL (> %.6f Mbps): %s\n', ...
    matThr, ternary(scan_el.anyMaterial,'SI','NO'));

%% ========================================================================
%  9. (C) BARRIDO DE ANCLAS DE ALPHA con las anclas de ELEVACION AJUSTADAS
%  ========================================================================
%  Es el barrido del punto 5 repetido con elMin/elMax = recorrido real. Responde a
%  la objecion exacta: el barrido de familia de E5 movia alpha_min/alpha_max con la
%  regla estrangulada, luego su conclusion ("ningun miembro de la familia gana")
%  podia ser un artefacto de las anclas de elevacion, no de la fisica.
fprintf('\n#######################################################################\n');
fprintf('#  (C) FAMILIA DE ANCLAS DE ALPHA con elevacion AJUSTADA (%.2f/%.2f deg)\n', ...
    elFitMin, elFitMax);
fprintf('#######################################################################\n');
cfg_c = cfg;  cfg_c.ffr.scheme = 'ffr';  cfg_c.ffr.Delta = DELTA;
cfg_c.ffr.adaptive = true;  cfg_c.ffr.adapt_dir = +1;
cfg_c.ffr.elMin_deg = elFitMin;  cfg_c.ffr.elMax_deg = elFitMax;
scan_alpha_anch = alpha_anchor_scan(cfg_c, ctx, sw, aLo, aHi);

fprintf('  MEJOR ganancia (R_agg, R_p5_borde): %+.6f Mbps  (anclas alpha %.2f-%.2f)\n', ...
    scan_alpha_anch.bestAgg, scan_alpha_anch.bestAgg_anchors(1), scan_alpha_anch.bestAgg_anchors(2));
fprintf('  MEJOR ganancia (R_p5,  R_p5_borde): %+.6f Mbps  (anclas alpha %.2f-%.2f)\n', ...
    scan_alpha_anch.bestP5, scan_alpha_anch.bestP5_anchors(1), scan_alpha_anch.bestP5_anchors(2));
fprintf('\n  COMPARACION con el barrido BASE (elevacion congelada en %g/%g):\n', EL_MIN_NOM, EL_MAX_NOM);
fprintf('    %-26s %14s %14s %14s\n', 'plano', 'base (25/90)', 'ajustadas', 'diferencia');
fprintf('    %-26s %+14.6f %+14.6f %+14.6f\n', '(R_agg,  R_p5_borde)', ...
    scan.bestAgg, scan_alpha_anch.bestAgg, scan_alpha_anch.bestAgg - scan.bestAgg);
fprintf('    %-26s %+14.6f %+14.6f %+14.6f\n', '(R_p5,   R_p5_borde)', ...
    scan.bestP5,  scan_alpha_anch.bestP5,  scan_alpha_anch.bestP5  - scan.bestP5);
fprintf('    umbral de materialidad: %.6f Mbps\n', matThr);
fprintf('    => con anclas ajustadas, la mejor de la familia %s el umbral material.\n', ...
    ternary(max(scan_alpha_anch.bestAgg, scan_alpha_anch.bestP5) > matThr, 'SUPERA', 'NO supera'));

%% ========================================================================
%  10. LINEALIDAD DE R_p5 DE BORDE EN alpha
%  ========================================================================
%  R_p5_borde(alpha) debe ser EXACTAMENTE lineal y pasar por el origen en alpha=1:
%  con sched='share' el ancho de borde es Be = (1-alpha)*B/Delta repartido entre
%  los usuarios de borde de la celda, y la SINR de borde NO depende de alpha (el
%  conjunto co-canal lo fija el color, no el reparto). Luego R = a*(1-alpha).
%  Que el residuo sea ~1e-12 es la prueba de que el "alpha equivalente" esta bien
%  definido: cualquier configuracion adaptativa tiene UN alpha fijo con su mismo
%  R_p5 de borde, y por eso la adaptacion solo puede DESPLAZARSE por la frontera.
fprintf('\n=============== LINEALIDAD DE R_p5 DE BORDE EN alpha ===============\n');
u   = 1 - alphas(:);
y   = sw.R_p5_edge(:);
aFit = (u.' * y) / (u.' * u);                  % minimos cuadrados sin termino independiente
res  = y - aFit*u;
linearity_residual = struct('a_fit', aFit, 'max_abs_res', max(abs(res)), ...
    'max_rel_res', max(abs(res))/max(abs(y)), 'A_EDGE_publicado', A_EDGE, ...
    'dif_a_vs_publicado', abs(aFit - A_EDGE), 'residuals', res.');
fprintf('  ajuste R = a*(1-alpha) sobre los %d puntos:  a = %.9f Mbps\n', numel(alphas), aFit);
fprintf('  residuo maximo |R - a*(1-alpha)| = %.3e Mbps (relativo %.3e)\n', ...
    linearity_residual.max_abs_res, linearity_residual.max_rel_res);
fprintf('  a publicado (usado en alpha_equiv) = %.6f  ->  |dif| = %.3e\n', ...
    A_EDGE, linearity_residual.dif_a_vs_publicado);
fprintf('  ==> R_p5 de BORDE es lineal en alpha a nivel de %s.\n', ...
    ternary(linearity_residual.max_abs_res < 1e-9, 'PRECISION DE MAQUINA', 'NO lineal: revisar'));

%% ========================================================================
%  11. (D) BARRIDO ESTATICO DE tau  (tau_q = 40 / 50 / 60)
%  ========================================================================
%  tau es ESTATICO en los tres bloques: lo unico que cambia es el CUANTIL fijo del
%  reparto centro/borde. No se adapta tau y no se toca ffr_policy.
fprintf('\n#######################################################################\n');
fprintf('#  (D) BARRIDO ESTATICO DE tau_q  (40 / 50 / 60), tau_mode = quantile\n');
fprintf('#######################################################################\n');
fprintf(['AVISO METODOLOGICO, y hay que respetarlo al leer la tabla: R_p5 de BORDE **NO es\n' ...
         'comparable entre tau distintos**, porque al mover el cuantil cambia la POBLACION\n' ...
         'de borde sobre la que se toma el percentil (y con ella el ancho por usuario). Las\n' ...
         'magnitudes que SI son comparables entre bloques son la GANANCIA SOBRE LA FRONTERA\n' ...
         '(cada una medida contra la frontera de SU propio barrido de alpha) y el CRITERIO\n' ...
         'Delta*SE_c/SE_b. No se cruzan valores absolutos entre bloques de tau.\n']);

tau_list  = [50 40 60];
tau_sweep = repmat(empty_tau_row(), 1, numel(tau_list));
tau_sweep(1) = pack_tau_row(TB50, ALPHA_MIN_NOM, ALPHA_MAX_NOM);
for it = 2:numel(tau_list)
    fprintf('\n#######################################################################\n');
    fprintf('#  BLOQUE tau_q = %d%%\n', tau_list(it));
    fprintf('#######################################################################\n');
    TBi = tau_block(cfg, ctx, maskV, alphas, DELTA, tau_list(it), true);
    tau_sweep(it) = pack_tau_row(TBi, ALPHA_MIN_NOM, ALPHA_MAX_NOM);
end

fprintf('\n================= RESUMEN DEL BARRIDO DE tau (magnitudes por bloque) =================\n');
fprintf('%6s %8s | %9s %10s | %9s %10s | %11s %10s %9s | %8s\n', ...
    'tau_q','alpha*','R5g_est','R5b_est','R5g_ada','R5b_ada','ganancia','umbral','material','critFFR');
fprintf('%s\n', repmat('-',1,116));
for it = 1:numel(tau_list)
    S = tau_sweep(it);
    fprintf('%5d%% %8.2f | %9.4f %10.4f | %9.4f %10.4f | %+11.6f %10.6f %9s | %8.6f\n', ...
        S.tau_q, S.alpha_star, S.R_p5_glob_static, S.R_p5_edge_static, ...
        S.R_p5_glob_adapt, S.R_p5_edge_adapt, S.gain_y, S.matThr, ...
        ternary(S.material,'SI','no'), S.ratio_all);
end
fprintf('%s\n', repmat('-',1,116));
fprintf(['  R5g/R5b = R_p5 global / R_p5 de BORDE. est = estatica alpha*; ada = adaptativa\n' ...
         '  nominal. "ganancia" = adaptativa sobre la frontera de compromiso de SU bloque.\n' ...
         '  RECORDATORIO: las columnas R5b NO se comparan entre filas (poblacion distinta).\n']);

fprintf('\n--- clasificacion centro/borde por bloque (control de que tau hace lo que dice) ---\n');
fprintf('%6s %14s %12s %12s %14s\n','tau_q','centerFrac','n_borde','n_cubiertas','pendiente a');
fprintf('%s\n', repmat('-',1,62));
for it = 1:numel(tau_list)
    S = tau_sweep(it);
    fprintf('%5d%% %14.6f %12d %12d %14.6f\n', ...
        S.tau_q, S.centerFrac, S.n_edge, S.n_cov, S.a_fit);
end
fprintf('%s\n', repmat('-',1,62));
fprintf(['  "pendiente a" = coeficiente de R_p5_borde = a*(1-alpha) en ESE bloque. Cambia\n' ...
         '  con tau (cambia la poblacion), y por eso el alpha equivalente solo se calcula\n' ...
         '  con la constante publicada %.6f en el bloque nominal tau_q = 50%%.\n'], A_EDGE);

fprintf('\n--- criterio Delta*SE_c/SE_b por tramo de elevacion y por bloque de tau ---\n');
fprintf('%6s | %-30s | %-30s | %-30s | %10s\n', 'tau_q', ...
    '25-35 deg (SEc/SEb/crit)','35-60 deg (SEc/SEb/crit)','60-90 deg (SEc/SEb/crit)','global');
fprintf('%s\n', repmat('-',1,125));
for it = 1:numel(tau_list)
    S = tau_sweep(it);
    fprintf('%5d%% |', S.tau_q);
    for e = 1:3
        if isnan(S.ratio_bin(e))
            fprintf(' %-30s |', '    -  (tramo vacio)');
        else
            fprintf(' %8.4f %8.4f %10.6f  |', S.bin_SE_c(e), S.bin_SE_b(e), S.ratio_bin(e));
        end
    end
    fprintf(' %10.6f\n', S.ratio_all);
end
fprintf('%s\n', repmat('-',1,125));

%% ========================================================================
%  12. PUERTA DE REGRESION (obligatoria ANTES de guardar)
%  ========================================================================
fprintf('\n#######################################################################\n');
fprintf('#  PUERTA DE REGRESION contra e5_adaptive_results.mat (tolerancia 1e-9)\n');
fprintf('#######################################################################\n');
TOLREG = 1e-9;

% LA REFERENCIA SE LEE DEL PROPIO .mat PUBLICADO, no de constantes transcritas a
% mano. Es una decision deliberada y es la que hace la puerta util: una constante
% copiada a 6 decimales no puede comprobarse a 1e-9 (su propio redondeo vale 5e-7),
% y ademas la transcripcion es justo donde se cuelan los errores. Contra el fichero
% la comparacion es a precision COMPLETA y no hay intermediario.
E5F = fullfile(ROOTDIR, 'e5_adaptive_results.mat');
if ~exist(E5F,'file')
    diary off;
    error('run_h2_anchors:sinRef', ...
        ['No se encuentra %s: sin el no hay puerta de regresion y este runner no ' ...
         'debe guardar nada.'], E5F);
end
E5R = load(E5F, 'KK','sw','GR','GRi','scan','alpha_star','ratio_all','ratio_bin', ...
                'alphas','DELTA');
fprintf('referencia: %s\n', E5F);

% El eje del barrido y el Delta tienen que ser los mismos o la comparacion no
% significa nada.
if ~isequal(E5R.alphas(:).', alphas(:).') || E5R.DELTA ~= DELTA
    diary off;
    error('run_h2_anchors:ejeRef', ...
        'El barrido de alpha o Delta del .mat de referencia no coinciden con los de este runner.');
end

% matThr NO se guarda en e5_adaptive_results.mat: se RECONSTRUYE con la misma
% formula a partir de sw y GR DEL FICHERO, asi que sigue anclado a la referencia.
resolY_ref = median(abs(diff(E5R.sw.R_p5_edge)));
matThr_ref = max(0.02*abs(E5R.GR.y_frontier), 0.5*resolY_ref);

REG = struct('name',{{}}, 'ref',[], 'got',[], 'dif',[], 'ok',[]);
% --- las CUATRO configuraciones originales, en sus tres KPIs de cabecera ---
for c = 1:4
    REG = reg_add(REG, sprintf('KK{%d}.all.R_p5_Mbps',  c), E5R.KK{c}.all.R_p5_Mbps,  KK{c}.all.R_p5_Mbps);
    REG = reg_add(REG, sprintf('KK{%d}.edge.R_p5_Mbps', c), E5R.KK{c}.edge.R_p5_Mbps, KK{c}.edge.R_p5_Mbps);
    REG = reg_add(REG, sprintf('KK{%d}.all.R_agg_Mbps', c), E5R.KK{c}.all.R_agg_Mbps, KK{c}.all.R_agg_Mbps);
end
% --- clasificacion centro/borde (si cambia, nada es comparable) ---
REG = reg_add(REG, 'KK{2}.centerFrac', E5R.KK{2}.centerFrac, KK{iST}.centerFrac);
REG = reg_add(REG, 'KK{2}.edge.n',     E5R.KK{2}.edge.n,     KK{iST}.edge.n);
REG = reg_add(REG, 'KK{2}.all.n',      E5R.KK{2}.all.n,      KK{iST}.all.n);
% --- el barrido de alpha COMPLETO, campo a campo (max|dif| sobre los 19 puntos) ---
swf = {'R_p5','R_agg','R_p5_edge','R_p5_V','R_agg_V','R_p5_edge_V'};
for i = 1:numel(swf)
    REG = reg_add_max(REG, sprintf('sw.%s (19 puntos, max|dif|)', swf{i}), ...
                      max(abs(sw.(swf{i}) - E5R.sw.(swf{i}))));
end
% --- escalares de decision ---
REG = reg_add(REG, 'alpha_star',              E5R.alpha_star, alpha_star);
REG = reg_add(REG, 'GR.gain_y',               E5R.GR.gain_y,  GR.gain_y);
REG = reg_add(REG, 'GRi.gain_y (invertida)',  E5R.GRi.gain_y, GRi.gain_y);
REG = reg_add(REG, 'umbral de materialidad',  matThr_ref,     matThr);
REG = reg_add(REG, 'scan.bestP5',             E5R.scan.bestP5, scan.bestP5);
REG = reg_add(REG, 'scan.bestP5 ancla inf.',  E5R.scan.bestP5_anchors(1), scan.bestP5_anchors(1));
REG = reg_add(REG, 'scan.bestP5 ancla sup.',  E5R.scan.bestP5_anchors(2), scan.bestP5_anchors(2));
REG = reg_add(REG, 'scan.bestAgg',            E5R.scan.bestAgg, scan.bestAgg);
REG = reg_add(REG, 'ratio_all',               E5R.ratio_all,  ratio_all);
% --- criterio por tramo, elemento a elemento (NaN == NaN cuenta como OK) ---
for e = 1:3
    if isnan(E5R.ratio_bin(e))
        REG = reg_add_max(REG, sprintf('ratio_bin(%d) [NaN vs NaN]', e), ...
                          double(~isnan(ratio_bin(e))));
    else
        REG = reg_add(REG, sprintf('ratio_bin(%d)', e), E5R.ratio_bin(e), ratio_bin(e));
    end
end

fprintf('\n%-40s %16s %16s %11s %s\n', 'magnitud','referencia','obtenido','|dif|','');
fprintf('%s\n', repmat('-',1,105));
nBad = 0;
for i = 1:numel(REG.name)
    d = abs(REG.got(i) - REG.ref(i));
    REG.dif(i) = d;
    REG.ok(i)  = d <= TOLREG;
    if ~REG.ok(i), nBad = nBad + 1; end
    fprintf('%-40s %16.9f %16.9f %11.2e %s\n', REG.name{i}, REG.ref(i), REG.got(i), d, ...
        ternary(REG.ok(i),'OK','*** FALLA ***'));
end
fprintf('%s\n', repmat('-',1,105));
REG.tol = TOLREG;  REG.nBad = nBad;
fprintf('%d de %d comprobaciones dentro de %.0e (max|dif| = %.2e).\n', ...
    sum(REG.ok), numel(REG.ok), TOLREG, max(REG.dif));

% --- SEGUNDA VUELTA: la tabla de cifras DECLARADAS en el encargo ---------
%  Se comprueban APARTE y a 1e-6 (su propio redondeo a 6 decimales no permite mas)
%  contra lo que dice el .mat. NO son la puerta: el .mat es la autoridad. Si una
%  discrepa, lo que hay es un error de TRANSCRIPCION en la cifra declarada, no una
%  regresion del codigo, y por eso se informa en vez de abortar. Es exactamente el
%  fallo contra el que existe esta segunda vuelta.
fprintf('\n--- Contraste de las cifras DECLARADAS en el encargo contra el .mat (1e-6) ---\n');
DECL = { ...
  'KK{1}.edge.R_p5_Mbps (reuse1)',      5.589765,  E5R.KK{1}.edge.R_p5_Mbps ; ...
  'KK{2}.edge.R_p5_Mbps (estatica)',    7.430042,  E5R.KK{2}.edge.R_p5_Mbps ; ...
  'KK{3}.edge.R_p5_Mbps (adaptativa)',  8.267286,  E5R.KK{3}.edge.R_p5_Mbps ; ...
  'KK{3}.all.R_p5_Mbps',                7.184932,  E5R.KK{3}.all.R_p5_Mbps  ; ...
  'alpha_star',                         0.55,      E5R.alpha_star           ; ...
  'GR.gain_y',                         -0.062665,  E5R.GR.gain_y            ; ...
  'GRi.gain_y (invertida)',            -0.213859,  E5R.GRi.gain_y           ; ...
  'umbral de materialidad',             0.412780,  matThr_ref               ; ...
  'sw.R_p5_edge(1)',                   15.685645,  E5R.sw.R_p5_edge(1)      ; ...
  'sw.R_p5_edge(19)',                   0.825560,  E5R.sw.R_p5_edge(end)    ; ...
  'scan.bestP5',                        0.053320,  E5R.scan.bestP5          ; ...
  'scan.bestP5 ancla inf.',             0.35,      E5R.scan.bestP5_anchors(1) ; ...
  'scan.bestP5 ancla sup.',             0.65,      E5R.scan.bestP5_anchors(2) ; ...
  'ratio_all',                          0.855612,  E5R.ratio_all            ; ...
  'ratio_bin(2)',                       0.777904,  E5R.ratio_bin(2)         ; ...
  'ratio_bin(3)',                       0.862729,  E5R.ratio_bin(3)         ; ...
  'centerFrac (tau=50%)',               0.500904,  E5R.KK{2}.centerFrac     ; ...
  'n_borde (tau=50%)',              33396,         E5R.KK{2}.edge.n         };
fprintf('%-38s %14s %18s %11s %s\n','cifra declarada','declarada','en el .mat','|dif|','');
fprintf('%s\n', repmat('-',1,100));
declBad = {};
for i = 1:size(DECL,1)
    d = abs(DECL{i,2} - DECL{i,3});
    ok = d <= 1e-6;
    if ~ok, declBad{end+1} = DECL{i,1}; end %#ok<SAGROW>
    fprintf('%-38s %14.6f %18.9f %11.2e %s\n', DECL{i,1}, DECL{i,2}, DECL{i,3}, d, ...
        ternary(ok,'OK','*** DISCREPA ***'));
end
fprintf('%s\n', repmat('-',1,100));
REG.declaradas_discrepantes = declBad;
if isempty(declBad)
    fprintf('Las %d cifras declaradas coinciden con el .mat.\n', size(DECL,1));
else
    fprintf(['AVISO: %d cifra(s) declarada(s) NO coinciden con e5_adaptive_results.mat:\n' ...
             '  %s\n' ...
             'NO es una regresion: este runner reproduce el .mat a %.0e en las %d\n' ...
             'comprobaciones de la puerta. Lo que discrepa es la cifra transcrita, y por eso\n' ...
             'la puerta se ancla al FICHERO y no a la transcripcion. Revisar el valor citado.\n'], ...
        numel(declBad), strjoin(declBad, ' | '), TOLREG, numel(REG.ok));
end

% --- Consistencia interna del criterio por tramo (identidad algebraica, 1e-12) ---
fprintf('\nConsistencia interna del criterio por tramo (bloque nominal, 1e-12):\n');
S1 = tau_sweep(1);
maxIntern = 0;
for e = 1:3
    if S1.bin_n(e) == 0
        fprintf('  tramo %d: VACIO (n = 0) -> ratio_bin = NaN\n', e);
        continue;
    end
    rr = DELTA * S1.bin_SE_c(e) / max(S1.bin_SE_b(e), eps);
    dI = abs(rr - S1.ratio_bin(e));
    maxIntern = max(maxIntern, dI);
    fprintf('  tramo %d: SE_c = %.6f | SE_b = %.6f | D*SEc/SEb = %.9f | |dif| = %.2e %s\n', ...
        e, S1.bin_SE_c(e), S1.bin_SE_b(e), rr, dI, ternary(dI<=1e-12,'OK','*** FALLA ***'));
end
if maxIntern > 1e-12
    nBad = nBad + 1;
    fprintf('  ==> DELTA*SE_c/SE_b no reproduce ratio_bin: los dos factores guardados no valen.\n');
end
REG.max_consistencia_bin = maxIntern;

% Linealidad de R_p5 de borde (control declarado, no de regresion).
fprintf('\nControl de linealidad: residuo maximo = %.3e Mbps (esperado ~1e-12) %s\n', ...
    linearity_residual.max_abs_res, ...
    ternary(linearity_residual.max_abs_res < 1e-9,'OK','*** REVISAR ***'));

fprintf('\n-----------------------------------------------------------------------\n');
if nBad > 0
    fprintf('PUERTA DE REGRESION FALLIDA: %d comprobacion(es) fuera de tolerancia.\n', nBad);
    fprintf('NO se guarda h2_anchors.mat.\n');
    diary off;
    error('run_h2_anchors:regresion', ...
        ['%d comprobacion(es) no reproducen e5_adaptive_results.mat dentro de %.0e. ' ...
         'Ver el log %s.'], nBad, TOLREG, LOGF);
end
fprintf('PUERTA DE REGRESION SUPERADA: E5 se reproduce y los bloques nuevos son comparables.\n');
fprintf('-----------------------------------------------------------------------\n');

%% ========================================================================
%  13. VEREDICTO
%  ========================================================================
fprintf('\n============================== VEREDICTO ==============================\n');
fprintf('(1) La regla NOMINAL esta ESTRANGULADA: usa el %.1f%% de su rango de alpha\n', ...
    100*anchors_used.frac_range_nom);
fprintf('    (alpha recorre [%.4f, %.4f] de un nominal [%.2f, %.2f], porque la elevacion\n', ...
    anchors_used.alpha_range_nom(1), anchors_used.alpha_range_nom(2), ALPHA_MIN_NOM, ALPHA_MAX_NOM);
fprintf('    media solo recorre [%.2f, %.2f] deg frente a las anclas %g/%g).\n', ...
    elFitMin, elFitMax, EL_MIN_NOM, EL_MAX_NOM);

fprintf('\n(2) Con las anclas AJUSTADAS al recorrido real (100%% del rango de alpha):\n');
fprintf('    ganancia sobre la frontera %+.6f Mbps (+1) y %+.6f Mbps (-1),\n', ...
    GR_anch.gain_y, GRi_anch.gain_y);
fprintf('    frente a %+.6f Mbps de la nominal. Umbral material: %.6f Mbps.\n', ...
    GR.gain_y, matThr);
fprintf('    => %s\n', ternary(GR_anch.gain_y > matThr, ...
    'AJUSTAR LAS ANCLAS SI DA MARGEN MATERIAL: el resultado de E5 estaba condicionado.', ...
    'NI ASI supera el umbral: el resultado negativo de H2 NO se debia a las anclas.'));

fprintf('\n(3) Barrido 2D de anclas de elevacion (16 combinaciones): mejor %+.6f Mbps\n', ...
    scan_el.bestP5);
fprintf('    en (elMin=%g, elMax=%g). Alguna material: %s.\n', ...
    scan_el.bestP5_anchors(1), scan_el.bestP5_anchors(2), ternary(scan_el.anyMaterial,'SI','NO'));

fprintf('\n(4) Familia de anclas de alpha con elevacion ajustada: mejor %+.6f Mbps\n', ...
    scan_alpha_anch.bestP5);
fprintf('    (base con elevacion congelada: %+.6f). %s\n', scan.bestP5, ...
    ternary(max(scan_alpha_anch.bestP5, scan_alpha_anch.bestAgg) > matThr, ...
    'SUPERA el umbral material.', 'Sigue SIN superar el umbral material.'));

fprintf('\n(5) Barrido estatico de tau (40/50/60): ganancia sobre la frontera de cada bloque\n');
fprintf('    %s Mbps; criterio Delta*SE_c/SE_b %s.\n', ...
    mat2str(round([tau_sweep.gain_y],6)), mat2str(round([tau_sweep.ratio_all],6)));
fprintf('    => H2 %s con ningun tau del barrido.\n', ...
    ternary(any([tau_sweep.material]), 'SE SOPORTA', 'NO se soporta'));

%% ========================================================================
%  14. GUARDADO
%  ========================================================================
timing = struct('warmup_p618_s', tWarm, 'experiment_s', toc(tExp), 'total_s', toc(tRun));
fprintf('\n[tiempos] experimento %.1f s (+ %.1f s de warmup P.618) | total %.1f s\n', ...
    timing.experiment_s, timing.warmup_p618_s, timing.total_s);

OUTF = fullfile(RESDIR, 'h2_anchors.mat');
save(OUTF, ...
     'cfg','BL','alphas','DELTA','elMin_valid','A_EDGE', ...
     'sw','GR','GRv','GRi','FR','FRv','FRi','scan','alpha_star','matThr', ...
     'ratio_bin','ratio_all','degenOK', ...
     'GR_anch','GRi_anch','FR_anch','FRi_anch', ...
     'scan_el','scan_alpha_anch','anchors_used','alpha_stats','linearity_residual', ...
     'tau_sweep','tau_list','REG','timing','MEM','BINFO');
fprintf('\nResultados guardados en %s\n', OUTF);
fprintf('Log en %s\n', LOGF);
fprintf('\nFIN.\n');
diary off;

%% ========================= funciones locales =========================

function [TB, cfg] = tau_block(cfg, ctx, maskV, alphas, DELTA, tau_q, verbose)
%TAU_BLOCK  Bloque COMPLETO de E5 para un cuantil tau fijo.
%   Barrido de alpha (frontera estatica) + alpha* + las 4 configuraciones +
%   frontier_gain en los dos planos + criterio por tramo de elevacion.
%
%   FUENTE UNICA: el bloque nominal (tau_q = 50) y los del barrido de tau salen de
%   AQUI, luego si el nominal reproduce E5 (puerta de regresion) los otros dos
%   recorren exactamente el mismo camino de codigo con otro cuantil.
%
%   El ORDEN de las mutaciones de cfg replica literalmente el de run_e5_adaptive:
%   el barrido deja cfg.ffr.alpha en su ultimo valor y los casos solo pisan alpha
%   y adapt_dir cuando el caso lo pide. Cambiarlo no deberia mover ningun numero,
%   pero replicarlo es lo que hace la reproduccion EXACTA y no aproximada.
%
%   Devuelve tambien la cfg mutada, para que el llamador continue en el mismo
%   estado en que run_e5_adaptive dejaba la suya.

cfg.ffr.tau_mode = 'quantile';
cfg.ffr.tau_q    = tau_q;

nA = numel(alphas);
sw = struct('alpha', alphas, 'Delta', DELTA, ...
            'R_p5', nan(1,nA), 'R_agg', nan(1,nA), 'R_p5_edge', nan(1,nA), ...
            'R_p5_V', nan(1,nA), 'R_agg_V', nan(1,nA), 'R_p5_edge_V', nan(1,nA));

cfg.ffr.scheme   = 'ffr';
cfg.ffr.Delta    = DELTA;
cfg.ffr.adaptive = false;

if verbose
    fprintf('\n--------- BARRIDO DE ALPHA FIJO (frontera estatica, Delta=%d, tau_q=%d%%) ---------\n', ...
        DELTA, tau_q);
    fprintf('%6s | %10s %12s %12s | %10s %12s %12s\n', ...
        'alpha','R_p5','R_p5_BORDE','R_agreg','R_p5(V)','R_p5_B(V)','R_agreg(V)');
    fprintf('%6s | %10s %12s %12s | %10s %12s %12s\n', ...
        '','[Mbps]','[Mbps]','[Gbps]','[Mbps]','[Mbps]','[Gbps]');
end
for ia = 1:nA
    cfg.ffr.alpha = alphas(ia);
    A = ffr_allocate(cfg, ctx, false);
    F = compute_sinr_ffr(cfg, ctx, A);
    K = compute_kpis(cfg, F, A);
    Kv= compute_kpis(cfg, F, A, maskV);
    sw.R_p5(ia)        = K.all.R_p5_Mbps;
    sw.R_agg(ia)       = K.all.R_agg_Mbps;
    sw.R_p5_edge(ia)   = K.edge.R_p5_Mbps;
    sw.R_p5_V(ia)      = Kv.all.R_p5_Mbps;
    sw.R_agg_V(ia)     = Kv.all.R_agg_Mbps;
    sw.R_p5_edge_V(ia) = Kv.edge.R_p5_Mbps;
    if verbose
        fprintf('%6.2f | %10.1f %12.1f %12.2f | %10.1f %12.1f %12.2f\n', ...
            alphas(ia), K.all.R_p5_Mbps, K.edge.R_p5_Mbps, K.all.R_agg_Mbps/1e3, ...
            Kv.all.R_p5_Mbps, Kv.edge.R_p5_Mbps, Kv.all.R_agg_Mbps/1e3);
    end
end

% alpha* = el de MEJOR COMPROMISO (max R_p5 global), mismo criterio que E5.
[~, iBest]  = max(sw.R_p5);
[~, iBestV] = max(sw.R_p5_V);
alpha_star  = alphas(iBest);
alpha_starV = alphas(iBestV);
if verbose
    fprintf('\n[estatica] alpha* (max R_p5 global, ventana completa) = %.2f\n', alpha_star);
end

% Las CUATRO configuraciones (fisica identica, solo cambia cfg.ffr)
cfg.ffr.scheme = 'ffr';  cfg.ffr.Delta = DELTA;
cases = {};
cases{end+1} = struct('name','reuse1',        'scheme','reuse1','Delta',1, ...
                      'adaptive',false,'alpha',NaN,       'dir',0);
cases{end+1} = struct('name','ffr-estatica',  'scheme','ffr',   'Delta',DELTA, ...
                      'adaptive',false,'alpha',alpha_star,'dir',0);
cases{end+1} = struct('name','ffr-adaptativa','scheme','ffr',   'Delta',DELTA, ...
                      'adaptive',true, 'alpha',NaN,       'dir',+1);
cases{end+1} = struct('name','ffr-adapt(inv)','scheme','ffr',   'Delta',DELTA, ...
                      'adaptive',true, 'alpha',NaN,       'dir',-1);
nC = numel(cases);
KK = cell(1,nC); KV = cell(1,nC); AA = cell(1,nC); FF = cell(1,nC);
for c = 1:nC
    cc = cases{c};
    cfg.ffr.scheme   = cc.scheme;
    cfg.ffr.Delta    = cc.Delta;
    cfg.ffr.adaptive = cc.adaptive;
    if ~isnan(cc.alpha), cfg.ffr.alpha = cc.alpha; end
    if cc.dir ~= 0,      cfg.ffr.adapt_dir = cc.dir; end
    AA{c} = ffr_allocate(cfg, ctx, false);
    FF{c} = compute_sinr_ffr(cfg, ctx, AA{c});
    KK{c} = compute_kpis(cfg, FF{c}, AA{c});
    KV{c} = compute_kpis(cfg, FF{c}, AA{c}, maskV);
end
iR1 = 1;  iST = 2;  iAD = 3;  iIN = 4;

if verbose
    printTable(cases, KK);
end

% Fronteras en los dos planos
FR  = frontier_gain(sw.R_agg/1e3,   sw.R_p5_edge,   KK{iAD}.all.R_agg_Mbps/1e3, KK{iAD}.edge.R_p5_Mbps);
FRv = frontier_gain(sw.R_agg_V/1e3, sw.R_p5_edge_V, KV{iAD}.all.R_agg_Mbps/1e3, KV{iAD}.edge.R_p5_Mbps);
FRi = frontier_gain(sw.R_agg/1e3,   sw.R_p5_edge,   KK{iIN}.all.R_agg_Mbps/1e3, KK{iIN}.edge.R_p5_Mbps);
FR.nPareto   = numel(FR.fx);
FR.degenPlane = (FR.nPareto <= 1);

GR  = frontier_gain(sw.R_p5,   sw.R_p5_edge,   KK{iAD}.all.R_p5_Mbps, KK{iAD}.edge.R_p5_Mbps);
GRv = frontier_gain(sw.R_p5_V, sw.R_p5_edge_V, KV{iAD}.all.R_p5_Mbps, KV{iAD}.edge.R_p5_Mbps);
GRi = frontier_gain(sw.R_p5,   sw.R_p5_edge,   KK{iIN}.all.R_p5_Mbps, KK{iIN}.edge.R_p5_Mbps);

resolY = median(abs(diff(sw.R_p5_edge)));
matThr = max(0.02*abs(GR.y_frontier), 0.5*resolY);
GR.material = GR.above && (GR.gain_y > matThr);

if verbose
    fprintf('\n--- frontera de compromiso (R_p5 global, R_p5 borde), tau_q=%d%% ---\n', tau_q);
    fprintf('  puntos NO dominados: %d de %d | plano (R_agg,R_borde) degenerado: %s (%d no dominados)\n', ...
        numel(GR.fx), nA, ternary(FR.degenPlane,'SI','no'), FR.nPareto);
    fprintf('  umbral de MATERIALIDAD: %.6f Mbps (2%% de la frontera / medio paso %.6f)\n', ...
        matThr, resolY);
    fprintf('  %-30s %8s %10s %14s %10s\n','','R_p5','R_p5_borde','frontera_borde','ganancia');
    report_frontier2('ADAPTATIVA  | ventana completa', GR);
    report_frontier2('ADAPTATIVA  | elev >= 45 deg  ', GRv);
    report_frontier2('adapt (inv) | ventana completa', GRi);
    fprintf('  correlacion(alpha(t), elevacion(t)) = %+.3f\n', ...
        corr_nan(AA{iAD}.alpha, AA{iAD}.POL.elevMean_deg));
end

% Criterio Delta*SE_c/SE_b por tramo de elevacion, con LOS DOS FACTORES aparte.
% Salen del MISMO compute_kpis con la mascara del tramo que ya calculaba E5; lo
% unico nuevo es que se guardan SE_c y SE_b y no solo su cociente.
edgesEl0 = [25 35 60 90];
ratio_bin = nan(1,3);  bin_SE_c = nan(1,3);  bin_SE_b = nan(1,3);  bin_n = nan(1,3);
if verbose
    fprintf('\n--- criterio Delta*SE_centro/SE_borde por tramo (tau_q=%d%%) ---\n', tau_q);
    fprintf('%-14s %10s %12s %12s %14s %8s\n', ...
        'tramo elev.','muestras','SE_centro','SE_borde','D*SEc/SEb','renta?');
end
for e = 1:3
    mk = ctx.cov & ctx.elev_deg >= edgesEl0(e) & ctx.elev_deg < edgesEl0(e+1);
    bin_n(e) = sum(mk(:));
    Kb = compute_kpis(cfg, FF{iST}, AA{iST}, mk);
    bin_SE_c(e) = Kb.center.SE_mean_bpsHz;
    bin_SE_b(e) = Kb.edge.SE_mean_bpsHz;
    ratio_bin(e) = DELTA * bin_SE_c(e) / max(bin_SE_b(e), eps);
    if verbose
        if bin_n(e) == 0
            fprintf('%-14s %10d %12s %12s %14s %8s\n', ...
                sprintf('%d-%d deg',edgesEl0(e),edgesEl0(e+1)), 0, 'n/d','n/d','NaN','n/d');
        else
            fprintf('%-14s %10d %12.6f %12.6f %14.6f %8s\n', ...
                sprintf('%d-%d deg',edgesEl0(e),edgesEl0(e+1)), bin_n(e), ...
                bin_SE_c(e), bin_SE_b(e), ratio_bin(e), ternary(ratio_bin(e)>1,'SI','NO'));
        end
    end
end
SE_c_all  = KK{iST}.center.SE_mean_bpsHz;
SE_b_all  = KK{iST}.edge.SE_mean_bpsHz;
ratio_all = DELTA * SE_c_all / max(SE_b_all, eps);
if verbose
    fprintf('%-14s %10d %12.6f %12.6f %14.6f %8s\n', 'TODO', KK{iST}.all.n, ...
        SE_c_all, SE_b_all, ratio_all, ternary(ratio_all>1,'SI','NO'));
end

% Pendiente de R_p5_borde = a*(1-alpha) EN ESTE BLOQUE (cambia con tau).
uu = 1 - alphas(:);
a_fit = (uu.' * sw.R_p5_edge(:)) / (uu.' * uu);

TB = struct('tau_q',tau_q, 'sw',sw, 'alphas',alphas, 'alpha_star',alpha_star, ...
    'alpha_starV',alpha_starV, 'cases',{cases}, 'KK',{KK}, 'KV',{KV}, ...
    'AA',{AA}, 'FF',{FF}, 'FR',FR, 'FRv',FRv, 'FRi',FRi, ...
    'GR',GR, 'GRv',GRv, 'GRi',GRi, 'matThr',matThr, 'resolY',resolY, ...
    'ratio_bin',ratio_bin, 'bin_SE_c',bin_SE_c, 'bin_SE_b',bin_SE_b, 'bin_n',bin_n, ...
    'ratio_all',ratio_all, 'SE_c_all',SE_c_all, 'SE_b_all',SE_b_all, ...
    'centerFrac',KK{iST}.centerFrac, 'n_edge',KK{iST}.edge.n, 'n_cov',KK{iST}.all.n, ...
    'a_fit',a_fit);
end

% -------------------------------------------------------------------------
function R = empty_tau_row()
%EMPTY_TAU_ROW  Plantilla de una fila del barrido de tau (para preasignar).
R = struct('tau_q',NaN,'alpha_star',NaN,'alphas',[],'sw',[], ...
    'R_p5_glob_static',NaN,'R_p5_edge_static',NaN,'R_agg_static',NaN, ...
    'R_p5_glob_adapt',NaN,'R_p5_edge_adapt',NaN,'R_agg_adapt',NaN, ...
    'R_p5_edge_reuse1',NaN,'R_p5_glob_reuse1',NaN, ...
    'gain_y',NaN,'gain_y_inv',NaN,'matThr',NaN,'material',false, ...
    'alpha_range_adapt',[NaN NaN],'frac_nominal_range',NaN, ...
    'ratio_bin',[NaN NaN NaN],'bin_SE_c',[NaN NaN NaN],'bin_SE_b',[NaN NaN NaN], ...
    'bin_n',[NaN NaN NaN],'ratio_all',NaN,'SE_c_all',NaN,'SE_b_all',NaN, ...
    'centerFrac',NaN,'n_edge',NaN,'n_cov',NaN,'a_fit',NaN,'GR',[]);
end

% -------------------------------------------------------------------------
function R = pack_tau_row(TB, aMinNom, aMaxNom)
%PACK_TAU_ROW  Reduce un bloque de tau a los escalares que se guardan.
%   Se descartan a proposito KK/KV/AA/FF: llevan las CDF completas (una por
%   muestra) y multiplicarian por 100 el tamano del .mat sin aportar nada que no
%   este ya resumido aqui.
R = empty_tau_row();
R.tau_q      = TB.tau_q;
R.alpha_star = TB.alpha_star;
R.alphas     = TB.alphas;
R.sw         = TB.sw;
R.R_p5_glob_reuse1 = TB.KK{1}.all.R_p5_Mbps;
R.R_p5_edge_reuse1 = TB.KK{1}.edge.R_p5_Mbps;
R.R_p5_glob_static = TB.KK{2}.all.R_p5_Mbps;
R.R_p5_edge_static = TB.KK{2}.edge.R_p5_Mbps;
R.R_agg_static     = TB.KK{2}.all.R_agg_Mbps;
R.R_p5_glob_adapt  = TB.KK{3}.all.R_p5_Mbps;
R.R_p5_edge_adapt  = TB.KK{3}.edge.R_p5_Mbps;
R.R_agg_adapt      = TB.KK{3}.all.R_agg_Mbps;
R.gain_y     = TB.GR.gain_y;
R.gain_y_inv = TB.GRi.gain_y;
R.matThr     = TB.matThr;
R.material   = TB.GR.above && (TB.GR.gain_y > TB.matThr);
av = TB.AA{3}.alpha;
R.alpha_range_adapt  = [min(av) max(av)];
R.frac_nominal_range = (max(av)-min(av)) / (aMaxNom - aMinNom);
R.ratio_bin = TB.ratio_bin;  R.bin_SE_c = TB.bin_SE_c;  R.bin_SE_b = TB.bin_SE_b;
R.bin_n     = TB.bin_n;      R.ratio_all = TB.ratio_all;
R.SE_c_all  = TB.SE_c_all;   R.SE_b_all  = TB.SE_b_all;
R.centerFrac = TB.centerFrac;  R.n_edge = TB.n_edge;  R.n_cov = TB.n_cov;
R.a_fit = TB.a_fit;
R.GR    = TB.GR;
end

% -------------------------------------------------------------------------
function S = alpha_anchor_scan(cfgS, ctx, sw, aLo, aHi)
%ALPHA_ANCHOR_SCAN  Barrido de la familia (alpha_min, alpha_max) sobre la frontera.
%   cfgS debe traer ya scheme/Delta/adaptive/adapt_dir y las anclas de ELEVACION
%   que se quieran usar: es lo unico que distingue el barrido base (25/90) del
%   barrido con anclas ajustadas. Fuente unica para los dos.
S = struct('aLo',aLo,'aHi',aHi, ...
           'gain_agg', -inf(numel(aLo),numel(aHi)), ...
           'gain_p5',  -inf(numel(aLo),numel(aHi)), ...
           'elMin_deg', cfgS.ffr.elMin_deg, 'elMax_deg', cfgS.ffr.elMax_deg);
for i = 1:numel(aLo)
    for j = 1:numel(aHi)
        if aHi(j) < aLo(i), continue; end
        cfgS.ffr.alpha_min = aLo(i);
        cfgS.ffr.alpha_max = aHi(j);
        As = ffr_allocate(cfgS, ctx, false);
        Fs = compute_sinr_ffr(cfgS, ctx, As);
        Ks = compute_kpis(cfgS, Fs, As);
        g1 = frontier_gain(sw.R_agg/1e3, sw.R_p5_edge, Ks.all.R_agg_Mbps/1e3, Ks.edge.R_p5_Mbps);
        g2 = frontier_gain(sw.R_p5,      sw.R_p5_edge, Ks.all.R_p5_Mbps,      Ks.edge.R_p5_Mbps);
        S.gain_agg(i,j) = g1.gain_y;
        S.gain_p5(i,j)  = g2.gain_y;
    end
end
[S.bestAgg, kA] = max(S.gain_agg(:));  [iA,jA] = ind2sub(size(S.gain_agg), kA);
[S.bestP5,  kP] = max(S.gain_p5(:));   [iP,jP] = ind2sub(size(S.gain_p5),  kP);
S.bestAgg_anchors = [aLo(iA) aHi(jA)];
S.bestP5_anchors  = [aLo(iP) aHi(jP)];
end

% -------------------------------------------------------------------------
function REG = reg_add(REG, name, ref, got)
%REG_ADD  Anade una comparacion escalar a la puerta de regresion.
REG.name{end+1} = name;
REG.ref(end+1)  = ref;
REG.got(end+1)  = got;
REG.dif(end+1)  = NaN;
REG.ok(end+1)   = false;
end

% -------------------------------------------------------------------------
function REG = reg_add_max(REG, name, maxdif)
%REG_ADD_MAX  Anade una comparacion VECTORIAL ya reducida a su max|dif|.
%   Se registra como referencia 0 frente a un "obtenido" que es la propia
%   desviacion, de modo que entra en la misma tabla y en el mismo umbral.
REG.name{end+1} = name;
REG.ref(end+1)  = 0;
REG.got(end+1)  = maxdif;
REG.dif(end+1)  = NaN;
REG.ok(end+1)   = false;
end

% -------------------------------------------------------------------------
function printTable(cases, KK)
fprintf('\n%-16s %8s %9s %9s %11s %9s %7s %7s\n', ...
    'configuracion','alpha','SINR_p5','R_p5','R_p5_BORDE','R_agreg','P(cob)','Jain');
fprintf('%-16s %8s %9s %9s %11s %9s %7s %7s\n', ...
    '','','[dB]','[Mbps]','[Mbps]','[Gbps]','','');
fprintf('%s\n', repmat('-',1,80));
for c = 1:numel(cases)
    K = KK{c};
    if strcmpi(cases{c}.scheme,'reuse1'), as = '   -   ';
    elseif cases{c}.adaptive,             as = ' adapt ';
    else,                                 as = sprintf('%7.2f', cases{c}.alpha);
    end
    fprintf('%-16s %8s %9.2f %9.4f %11.4f %9.4f %7.2f %7.3f\n', ...
        cases{c}.name, as, K.all.SINR_p5, K.all.R_p5_Mbps, K.edge.R_p5_Mbps, ...
        K.all.R_agg_Mbps/1e3, K.all.Pcov, K.all.jain);
end
fprintf('%s\n', repmat('-',1,80));
end

% -------------------------------------------------------------------------
% COPIA de la funcion local de run_e5_adaptive.m (no invocable desde fuera).
% Anclada por la puerta de regresion: GR.gain_y, GRi.gain_y y scan.bestP5 salen
% de aqui, luego si la copia divergiera el script abortaria.
function FR = frontier_gain(x, y, xq, yq)
%FRONTIER_GAIN  Posicion de (xq,yq) respecto a la FRONTERA DE PARETO de (x,y),
%   con el criterio "mas es mejor" en AMBOS ejes.
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
%INTERP_CLAMPED  Interpolacion lineal con extrapolacion PLANA en los extremos.
if numel(xx) < 2
    if isempty(xx), yq = NaN; else, yq = yy(1); end
    return;
end
if xq <= xx(1),   yq = yy(1);   return; end
if xq >= xx(end), yq = yy(end); return; end
yq = interp1(xx, yy, xq, 'linear');
end

% -------------------------------------------------------------------------
function report_frontier2(label, FR)
fprintf('  %-30s %8.4f %10.4f %14.4f %+10.6f  %s\n', ...
    label, FR.x_q, FR.y_q, FR.y_frontier, FR.gain_y, ternary(FR.above,'ENCIMA','debajo'));
end

% -------------------------------------------------------------------------
function r = corr_nan(a, b)
ok = isfinite(a) & isfinite(b);
if sum(ok) < 3 || std(a(ok)) == 0 || std(b(ok)) == 0, r = NaN; return; end
c = corrcoef(a(ok), b(ok));  r = c(1,2);
end

% -------------------------------------------------------------------------
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
