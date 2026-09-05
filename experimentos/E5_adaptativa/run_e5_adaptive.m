%% RUN_E5_ADAPTIVE  Experimento E5: FFR ADAPTATIVA (H2) frente a FFR estatica.
%  RUNNER DELGADO. Toda la fisica esta en el MOTOR y es EXACTAMENTE la misma que
%  usa run_ffr_demo (y que E0): aqui solo se varia cfg.
%
%     build_beam_layout -> ffr_context -> ffr_coloring -> ffr_policy
%                       -> ffr_allocate -> compute_sinr_ffr -> compute_kpis
%
%  HIPOTESIS H2: la FFR adaptativa recupera eficiencia espectral en geometrias
%  favorables SIN perder proteccion de borde. Es decir, su punto (R_agregado,
%  R_p5_borde) debe caer POR ENCIMA de la FRONTERA DE COMPROMISO que trazan todos
%  los alpha FIJOS posibles. Si cae SOBRE la frontera, la adaptacion no aporta.
%
%  REGLA ADAPTATIVA (implementada en ffr_policy, ver alli la justificacion):
%     alpha(t) = alpha_min + (alpha_max-alpha_min)*clip((elev(t)-elMin)/(elMax-elMin),0,1)
%  Como alpha es la fraccion del INTERIOR, alpha CRECIENTE con la elevacion
%  significa: baja elevacion -> banda al BORDE (protege); alta elevacion -> banda
%  al INTERIOR (libera agregado). Es el sentido MEDIDO en el barrido de alpha fijo.
%
%  Comparacion (misma poblacion de usuarios, mismos instantes, misma fisica):
%     reuse1 | ffr-estatica (mejor alpha fijo) | ffr-adaptativa | barrido completo
%
%  Uso: ejecutar en la raiz del proyecto.  Resultados en e5_adaptive_results.mat.

clear; clc; close all;

% CRONOMETRO DE EXTREMO A EXTREMO. El warmup de la
% tabla P.618 se mide APARTE y se descuenta: coste fijo de la configuracion
% atmosferica, cacheado en 'persistent' e independiente del experimento.
tRun   = tic;
FIGDIR = 'figs_e5';

%% 1. Perfil de estudio: IDENTICO al de run_ffr_demo (comparabilidad directa)
cfg = config_default();

% CONSTELACION REAL Starlink Shell-1 (Walker 53:1584/72/1). Antes T=66/P=6; el
% cambio va acompasado con run_ffr_demo para que E2 y E5 sigan siendo directamente
% comparables (que es la razon de que compartan perfil). Ver la nota extensa en
% run_ffr_demo.
cfg.constellations(1).T = 1584;
cfg.constellations(1).P = 72;

cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 3;

cfg.time.dt       = 60;            % s
cfg.time.duration = 7200;          % s (2 h: recorre elevaciones de 25 a ~81 deg)

cfg.ffr.tau_mode = 'quantile';     % reparto centro/borde fijo al 50%: asi la unica
cfg.ffr.tau_q    = 50;             % diferencia entre estatica y adaptativa es ALPHA

cfg.beams.nRings = 5;              % 91 celdas Earth-fixed (>=1.5 anillos de guarda)

% TROCEADO TEMPORAL: obligatorio a densidad real (7.53 GB de una pieza frente a
% 7.45 GB de RAM total). EXACTO, verificado en test_ctx_blocked (max|dif| = 0
% incluso contra el .mat publicado). Ver run_ffr_demo.
cfg.compute.timeBlock = 10;                                        % <-- DECISION

DELTA       = 3;                   % Delta de cabecera del experimento
elMin_valid = 45;                  % deg  submuestra de "geometria valida"
alphas      = 0.05:0.05:0.95;      % barrido fino: traza la FRONTERA de compromiso

fprintf('===== E5: FFR ADAPTATIVA (H2) =====\n');
fprintf('Ku DL: f=%.0f GHz | B=%.0f MHz | EIRPdens=%.1f dBW/MHz | HPBW=%.2f deg | G/T=%.1f dB/K\n', ...
    cfg.radio.freq_GHz, cfg.radio.B_MHz, cfg.radio.EIRPdensity_dBWMHz, ...
    cfg.radio.beamwidth3dB_deg, cfg.radio.GT_dBK);

% Warmup EXPLICITO de la tabla P.618 (coste fijo, fuera del cronometro). No cambia
% ningun numero: test_atm_cache verifica que la cache es identica bit a bit.
tW = tic;  atm_loss_dB(cfg, 45);  tWarm = toc(tW);
fprintf('[warmup] tabla P.618 construida en %.1f s (coste fijo, se descuenta)\n', tWarm);
tExp = tic;                        % <- cronometro del EXPERIMENTO, ya con cache

%% 2. Dimensionado de memoria + pipeline por BLOQUES TEMPORALES
%  Mismo camino que run_ffr_demo: build_ctx_blocked devuelve el ctx COMPLETO, luego
%  todo el analisis de E5 (frontera de Pareto, regla adaptativa, validaciones) es
%  identico al de antes. El troceado no cambia ningun numero (test_ctx_blocked).
MEM = estimate_sweep_memory(cfg, {struct('T', cfg.constellations(1).T)}, 1, ...
                            struct('verbose', true));

[ctx, sats, users, BL, BINFO] = build_ctx_blocked(cfg, true);
tvec = ctx.tvec;

fprintf('Rejilla local: %d usuarios (radio %.0f km, paso %.0f km) | Nt=%d | N=%d satelites\n', ...
    users.M, cfg.ground.radius_km, cfg.ground.step_km, ctx.Nt, sats.N);
fprintf('\nLayout: %d haces | s = %.4f deg = %.2f km\n', BL.nBeams, BL.spacing_deg, BL.spacing_km);

% SINR de referencia (reuso-1): la calcula UNA vez ffr_allocate y se cachea en ctx
% para que TODOS los esquemas clasifiquen centro/borde con la MISMA referencia.
cfg.ffr.scheme = 'reuse1';  cfg.ffr.Delta = 1;  cfg.ffr.adaptive = false;
A_r1 = ffr_allocate(cfg, ctx);
ctx.SINR_ref_dB = A_r1.SINR_ref_dB;
F_r1 = compute_sinr_ffr(cfg, ctx, A_r1);

maskV = ctx.cov & ctx.elev_deg >= elMin_valid;

%% 3. Barrido de ALPHA FIJO -> frontera de compromiso (referencia de la estatica)
nA = numel(alphas);
sw = struct('alpha', alphas, 'Delta', DELTA, ...
            'R_p5', nan(1,nA), 'R_agg', nan(1,nA), 'R_p5_edge', nan(1,nA), ...
            'R_p5_V', nan(1,nA), 'R_agg_V', nan(1,nA), 'R_p5_edge_V', nan(1,nA));
AA_sw = cell(1,nA);  FF_sw = cell(1,nA);

cfg.ffr.scheme   = 'ffr';
cfg.ffr.Delta    = DELTA;
cfg.ffr.adaptive = false;

fprintf('\n--------- BARRIDO DE ALPHA FIJO (frontera estatica, Delta=%d) ---------\n', DELTA);
fprintf('%6s | %10s %12s %12s | %10s %12s %12s\n', ...
    'alpha','R_p5','R_p5_BORDE','R_agreg','R_p5(V)','R_p5_B(V)','R_agreg(V)');
fprintf('%6s | %10s %12s %12s | %10s %12s %12s\n', ...
    '','[Mbps]','[Mbps]','[Gbps]','[Mbps]','[Mbps]','[Gbps]');
for ia = 1:nA
    cfg.ffr.alpha = alphas(ia);
    A = ffr_allocate(cfg, ctx, false);
    F = compute_sinr_ffr(cfg, ctx, A);
    K = compute_kpis(cfg, F, A);
    Kv= compute_kpis(cfg, F, A, maskV);
    AA_sw{ia} = A;  FF_sw{ia} = F;
    sw.R_p5(ia)        = K.all.R_p5_Mbps;
    sw.R_agg(ia)       = K.all.R_agg_Mbps;
    sw.R_p5_edge(ia)   = K.edge.R_p5_Mbps;
    sw.R_p5_V(ia)      = Kv.all.R_p5_Mbps;
    sw.R_agg_V(ia)     = Kv.all.R_agg_Mbps;
    sw.R_p5_edge_V(ia) = Kv.edge.R_p5_Mbps;
    fprintf('%6.2f | %10.1f %12.1f %12.2f | %10.1f %12.1f %12.2f\n', ...
        alphas(ia), K.all.R_p5_Mbps, K.edge.R_p5_Mbps, K.all.R_agg_Mbps/1e3, ...
        Kv.all.R_p5_Mbps, Kv.edge.R_p5_Mbps, Kv.all.R_agg_Mbps/1e3);
end

% alpha* de la estatica = el de MEJOR COMPROMISO (max R_p5 global), mismo criterio
% que en run_ffr_demo, para que la "ffr-estatica" de E5 sea la mejor estatica posible.
[~, iBest]  = max(sw.R_p5);
[~, iBestV] = max(sw.R_p5_V);
alpha_star  = alphas(iBest);
alpha_starV = alphas(iBestV);
fprintf('\n[estatica] alpha* (max R_p5 global, ventana completa)  = %.2f\n', alpha_star);
fprintf('[estatica] alpha* (max R_p5 global, elev >= %d deg)     = %.2f\n', elMin_valid, alpha_starV);

%% 4. Las CUATRO configuraciones (fisica identica, solo cambia cfg.ffr)
cfg.ffr.scheme = 'ffr';  cfg.ffr.Delta = DELTA;

cases = {};
cases{end+1} = struct('name','reuse1',        'scheme','reuse1','Delta',1, ...
                      'adaptive',false,'alpha',NaN,       'dir',0);
cases{end+1} = struct('name','ffr-estatica',  'scheme','ffr',   'Delta',DELTA, ...
                      'adaptive',false,'alpha',alpha_star,'dir',0);
cases{end+1} = struct('name','ffr-adaptativa','scheme','ffr',   'Delta',DELTA, ...
                      'adaptive',true, 'alpha',NaN,       'dir',+1);
% Contraste de SENTIDO (no es un competidor: es la prueba de que la ganancia viene
% de adaptar en la direccion CORRECTA y no del hecho de variar alpha por variar).
cases{end+1} = struct('name','ffr-adapt(inv)','scheme','ffr',   'Delta',DELTA, ...
                      'adaptive',true, 'alpha',NaN,       'dir',-1);
nC = numel(cases);
KK = cell(1,nC); KV = cell(1,nC); AA = cell(1,nC); FF = cell(1,nC);

fprintf('\n--------- EVALUACION POR CONFIGURACION ---------\n');
for c = 1:nC
    cc = cases{c};
    cfg.ffr.scheme   = cc.scheme;
    cfg.ffr.Delta    = cc.Delta;
    cfg.ffr.adaptive = cc.adaptive;
    if ~isnan(cc.alpha), cfg.ffr.alpha = cc.alpha; end
    if cc.dir ~= 0,      cfg.ffr.adapt_dir = cc.dir; end

    fprintf('\n>> %s\n', cc.name);
    AA{c} = ffr_allocate(cfg, ctx);
    FF{c} = compute_sinr_ffr(cfg, ctx, AA{c});
    KK{c} = compute_kpis(cfg, FF{c}, AA{c});
    KV{c} = compute_kpis(cfg, FF{c}, AA{c}, maskV);
    if cc.adaptive
        P = AA{c}.POL;
        fprintf(['   [politica] adapt_dir=%+d | alpha en [%.3f, %.3f] (anclas %.1f-%.1f deg, ' ...
                 'rango pedido %.2f-%.2f) | elev media %.1f-%.1f deg\n'], ...
            P.adapt_dir, P.alpha_range(1), P.alpha_range(2), P.elMin_deg, P.elMax_deg, ...
            P.alpha_min, P.alpha_max, min(P.elevMean_deg), max(P.elevMean_deg));
    end
end

iR1 = 1;  iST = 2;  iAD = 3;  iIN = 4;

%% 5. (a) Tabla: R_p5 de borde, R_agregado y R_p5 global
fprintf('\n==================== E5 | VENTANA COMPLETA (%d muestras) ====================\n', ...
    sum(ctx.cov(:)));
printTable(cases, KK);
fprintf('\n============ E5 | GEOMETRIA VALIDA (elev >= %d deg, %d muestras) ============\n', ...
    elMin_valid, sum(maskV(:)));
printTable(cases, KV);

fprintf('\n[H2] R_p5 de BORDE: reuse1 %.1f -> estatica(a=%.2f) %.1f -> ADAPTATIVA %.1f Mbps\n', ...
    KK{iR1}.edge.R_p5_Mbps, alpha_star, KK{iST}.edge.R_p5_Mbps, KK{iAD}.edge.R_p5_Mbps);
fprintf('[H2] R_agregado:    reuse1 %.2f -> estatica %.2f -> ADAPTATIVA %.2f Gbps\n', ...
    KK{iR1}.all.R_agg_Mbps/1e3, KK{iST}.all.R_agg_Mbps/1e3, KK{iAD}.all.R_agg_Mbps/1e3);

%% 6. (b) FRONTERA DE COMPROMISO: (R_agregado, R_p5_borde)
%  Prueba de H2: el punto adaptativo debe quedar POR ENCIMA de la frontera que
%  trazan los alpha fijos. "Por encima" = a igual R_agregado, mas R_p5 de borde.
fprintf('\n============== FRONTERA DE COMPROMISO (R_agregado, R_p5_borde) ==============\n');
FR   = frontier_gain(sw.R_agg/1e3,   sw.R_p5_edge,   KK{iAD}.all.R_agg_Mbps/1e3, KK{iAD}.edge.R_p5_Mbps);
FRv  = frontier_gain(sw.R_agg_V/1e3, sw.R_p5_edge_V, KV{iAD}.all.R_agg_Mbps/1e3, KV{iAD}.edge.R_p5_Mbps);
FRi  = frontier_gain(sw.R_agg/1e3,   sw.R_p5_edge,   KK{iIN}.all.R_agg_Mbps/1e3, KK{iIN}.edge.R_p5_Mbps);

report_frontier('ADAPTATIVA  | ventana completa', FR);
report_frontier('ADAPTATIVA  | elev >= 45 deg  ', FRv);
report_frontier('adapt (inv) | ventana completa', FRi);

% -- DIAGNOSTICO CRITICO: existe siquiera un compromiso en este plano? --
%  Si R_agregado Y R_p5_borde son AMBOS monotonos decrecientes en alpha, entonces
%  alpha -> 0 DOMINA en las dos metricas: no hay curva de compromiso que superar,
%  la "frontera" degenera en UN SOLO punto (el alpha mas pequeno) y la pregunta
%  "cae la adaptativa por encima?" no tiene contenido. Hay que detectarlo y decirlo.
FR.nPareto  = numel(FR.fx);
FR.degenPlane = (FR.nPareto <= 1);
monoAgg  = all(diff(sw.R_agg)  < 0);
monoEdge = all(diff(sw.R_p5_edge) < 0);
fprintf(['\n  [diagnostico del plano] puntos NO dominados del barrido: %d de %d.\n' ...
         '  R_agregado monotono DECRECIENTE en alpha: %s | R_p5_borde monotono DECRECIENTE: %s\n'], ...
    FR.nPareto, nA, ternary(monoAgg,'SI','no'), ternary(monoEdge,'SI','no'));
if FR.degenPlane
    fprintf(['  ==> El plano (R_agregado, R_p5_borde) NO CONTIENE UN COMPROMISO en este\n' ...
             '      escenario: dar banda al INTERIOR empeora AMBAS metricas a la vez, luego\n' ...
             '      alpha->0 domina y la frontera se reduce a un punto. La prueba "por encima\n' ...
             '      de la frontera" es aqui VACUA; el compromiso real esta en el plano\n' ...
             '      (R_p5 GLOBAL, R_p5 BORDE), que es el que se evalua a continuacion.\n']);
end

%% 6b. FRONTERA REAL DE ESTE SISTEMA: (R_p5 GLOBAL, R_p5 BORDE)
%  Es el compromiso que de verdad existe aqui y el que fijo alpha* en la FFR
%  estatica: R_p5 GLOBAL tiene un MAXIMO interior en alpha (si alpha es muy bajo,
%  los usuarios de CENTRO se quedan sin banda y hunden el percentil 5 global),
%  mientras R_p5 de BORDE es monotono decreciente en alpha. Proteger el borde
%  cuesta, por tanto, percentil global: ESE es el compromiso que la FFR adaptativa
%  tiene que batir para sostener H2.
fprintf('\n========= FRONTERA DE COMPROMISO REAL (R_p5 global, R_p5 borde) =========\n');
GR  = frontier_gain(sw.R_p5,   sw.R_p5_edge,   KK{iAD}.all.R_p5_Mbps, KK{iAD}.edge.R_p5_Mbps);
GRv = frontier_gain(sw.R_p5_V, sw.R_p5_edge_V, KV{iAD}.all.R_p5_Mbps, KV{iAD}.edge.R_p5_Mbps);
GRi = frontier_gain(sw.R_p5,   sw.R_p5_edge,   KK{iIN}.all.R_p5_Mbps, KK{iIN}.edge.R_p5_Mbps);
fprintf('  puntos NO dominados del barrido: %d de %d (si es >1, el compromiso EXISTE)\n', ...
    numel(GR.fx), nA);

% MATERIALIDAD de la ganancia. "Estar por encima" con signo positivo NO basta: una
% ganancia de milesimas de Mbps es indistinguible de la resolucion del propio
% barrido y no soporta H2. Se exige superar a la vez un umbral RELATIVO (2% del
% valor de la frontera) y la mitad del paso en y del barrido de alpha.
resolY = median(abs(diff(sw.R_p5_edge)));
matThr = max(0.02*abs(GR.y_frontier), 0.5*resolY);
GR.material = GR.above && (GR.gain_y > matThr);
fprintf('  umbral de MATERIALIDAD: %.3f Mbps (2%% de la frontera / medio paso del barrido %.3f)\n', ...
    matThr, resolY);
fprintf('  %-30s %8s %10s %12s %10s\n','','R_p5','R_p5_borde','frontera_borde','ganancia');
report_frontier2('ADAPTATIVA  | ventana completa', GR);
report_frontier2('ADAPTATIVA  | elev >= 45 deg  ', GRv);
report_frontier2('adapt (inv) | ventana completa', GRi);

%% 6c. POR QUE: criterio analitico de rentabilidad de la sub-banda INTERIOR
%  Con la convencion 'share' el caudal agregado por celda es
%       R_celda(alpha) = Bc*SE_centro + Be*SE_borde,   Bc = alpha*B, Be = (1-alpha)*B/Delta
%  luego                d R_celda / d alpha  =  B * ( SE_centro - SE_borde/Delta ).
%  La sub-banda INTERIOR solo RENTA (subir alpha aumenta el agregado) si
%
%       SE_centro > SE_borde / Delta      <=>      Delta * SE_centro / SE_borde > 1
%
%  Esta es la CONDICION DE VIABILIDAD de la FFR adaptativa: H2 postula que en
%  geometria favorable el interior rinde y se le puede devolver banda. Si el
%  cociente es < 1 en TODO el rango de elevaciones simulado, no hay ninguna
%  geometria en la que subir alpha compense y la adaptacion NO puede ganar:
%  el resultado negativo seria ESTRUCTURAL, no un problema de ajuste de anclas.
fprintf('\n===== CONDICION DE RENTABILIDAD DEL INTERIOR: Delta*SE_centro/SE_borde > 1 =====\n');
fprintf('%-14s %12s %12s %14s %10s\n', ...
    'tramo elev.','SE_centro','SE_borde','D*SEc/SEb','renta?');
ratio_bin = nan(1,3);
edgesEl0 = [25 35 60 90];
for e = 1:3
    mk = ctx.cov & ctx.elev_deg >= edgesEl0(e) & ctx.elev_deg < edgesEl0(e+1);
    Kb = compute_kpis(cfg, FF{iST}, AA{iST}, mk);
    rr = DELTA * Kb.center.SE_mean_bpsHz / max(Kb.edge.SE_mean_bpsHz, eps);
    ratio_bin(e) = rr;
    fprintf('%-14s %12.4f %12.4f %14.3f %10s\n', ...
        sprintf('%d-%d deg',edgesEl0(e),edgesEl0(e+1)), ...
        Kb.center.SE_mean_bpsHz, Kb.edge.SE_mean_bpsHz, rr, ternary(rr>1,'SI','NO'));
end
Kfull = KK{iST};
ratio_all = DELTA * Kfull.center.SE_mean_bpsHz / max(Kfull.edge.SE_mean_bpsHz, eps);
fprintf('%-14s %12.4f %12.4f %14.3f %10s\n', 'TODO', ...
    Kfull.center.SE_mean_bpsHz, Kfull.edge.SE_mean_bpsHz, ratio_all, ternary(ratio_all>1,'SI','NO'));

%% 7. (c) alpha(t) y elevacion: la regla en accion
POLa = AA{iAD}.POL;
elev_t = POLa.elevMean_deg;
alp_t  = AA{iAD}.alpha;
fprintf('\n[regla] correlacion(alpha(t), elevacion(t)) = %+.3f  (esperado ~ +1 con adapt_dir=+1)\n', ...
    corr_nan(alp_t, elev_t));

%% 8. (d) Desglose por TRAMOS DE ELEVACION: donde aporta mas la adaptacion
edgesEl = [25 35 60 90];
nB      = numel(edgesEl)-1;
binName = cell(1,nB);
for e = 1:nB, binName{e} = sprintf('%d-%d deg', edgesEl(e), edgesEl(e+1)); end

bin = struct('name',{binName}, 'n',nan(1,nB), 'alphaMean',nan(1,nB), ...
             'R_p5_edge',nan(nC,nB), 'R_agg',nan(nC,nB), 'R_p5',nan(nC,nB), ...
             'SINR_p5_edge',nan(nC,nB));
fprintf('\n=============== DESGLOSE POR TRAMO DE ELEVACION DEL SERVIDOR ===============\n');
for e = 1:nB
    mk = ctx.cov & ctx.elev_deg >= edgesEl(e) & ctx.elev_deg < edgesEl(e+1);
    bin.n(e) = sum(mk(:));
    kk = any(mk,1);
    bin.alphaMean(e) = mean(alp_t(kk));
    for c = 1:nC
        Kb = compute_kpis(cfg, FF{c}, AA{c}, mk);
        bin.R_p5_edge(c,e)    = Kb.edge.R_p5_Mbps;
        bin.R_agg(c,e)        = Kb.all.R_agg_Mbps;
        bin.R_p5(c,e)         = Kb.all.R_p5_Mbps;
        bin.SINR_p5_edge(c,e) = Kb.edge.SINR_p5;
    end
end

fprintf('%-16s', 'tramo');            for e=1:nB, fprintf('%14s', binName{e}); end
fprintf('\n%-16s', 'muestras');       for e=1:nB, fprintf('%14d', bin.n(e)); end
fprintf('\n%-16s', 'alpha adapt.');   for e=1:nB, fprintf('%14.3f', bin.alphaMean(e)); end
fprintf('\n%s\n', repmat('-',1,16+14*nB));
fprintf('R_p5 de BORDE [Mbps]:\n');
for c = 1:nC
    fprintf('  %-14s', cases{c}.name);
    for e=1:nB, fprintf('%14.1f', bin.R_p5_edge(c,e)); end
    fprintf('\n');
end
fprintf('  %-14s', 'GANANCIA adap.');
for e=1:nB
    g = bin.R_p5_edge(iAD,e) - bin.R_p5_edge(iST,e);
    fprintf('%13.1f%%', 100*g/max(bin.R_p5_edge(iST,e),eps));
end
fprintf('\n%s\nR_agregado [Gbps]:\n', repmat('-',1,16+14*nB));
for c = 1:nC
    fprintf('  %-14s', cases{c}.name);
    for e=1:nB, fprintf('%14.2f', bin.R_agg(c,e)/1e3); end
    fprintf('\n');
end
fprintf('  %-14s', 'GANANCIA adap.');
for e=1:nB
    g = bin.R_agg(iAD,e) - bin.R_agg(iST,e);
    fprintf('%13.1f%%', 100*g/max(bin.R_agg(iST,e),eps));
end
fprintf('\n%s\n', repmat('-',1,16+14*nB));

%% 9. TAREA 3 | VALIDACION: caso DEGENERADO alpha_min == alpha_max
%  Si la regla adaptativa se colapsa a un alpha constante, sus KPIs deben coincidir
%  EXACTAMENTE con los de la estatica con ese mismo alpha. Es la prueba de que la
%  adaptacion no introduce ninguna ventaja espuria (un bug que la favorezca).
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
fprintf('  max|d alpha|      = %.3g\n', dAlpha);
fprintf('  max|d SINR|       = %.3g dB\n', dSINR);
fprintf('  max|d B_user|     = %.3g Hz\n', dBu);
fprintf('  d clasificacion   = %d muestras\n', dCls);
fprintf('  max|d KPI|        = %.3g Mbps\n', dKpi);
degenOK = (dAlpha==0) && (dSINR==0) && (dBu==0) && (dCls==0) && (dKpi==0);
if degenOK
    fprintf('  ==> OK: la adaptativa degenera EXACTAMENTE en la estatica. Sin ventaja espuria.\n');
else
    warning('run_e5_adaptive:degenerate', ...
        'El caso degenerado NO coincide con la estatica: hay un bug en la politica adaptativa.');
end

%% 9b. Es el resultado ESTRUCTURAL o solo unas anclas mal elegidas?
%  DIAGNOSTICO, no ajuste: se recorre la FAMILIA de anclas (alpha_min, alpha_max) y
%  se mide la MEJOR ganancia alcanzable sobre la frontera estatica. Si ni el mejor
%  miembro de la familia supera la frontera, el resultado negativo NO se debe a la
%  eleccion de anclas (0.3/0.7) sino a la fisica del escenario.
%  NOTA: el maximo de este barrido NO se reporta como resultado de H2 (seria
%  ajustar las anclas para ganar); solo se usa para acotar el margen de la familia.
fprintf('\n===== DIAGNOSTICO: lo arreglarian otras anclas? (barrido de la FAMILIA) =====\n');
aLo = 0.05:0.10:0.55;
aHi = 0.05:0.10:0.95;
scan = struct('aLo',aLo,'aHi',aHi, ...
              'gain_agg', -inf(numel(aLo),numel(aHi)), ...
              'gain_p5',  -inf(numel(aLo),numel(aHi)));
cfg_s = cfg;  cfg_s.ffr.scheme = 'ffr';  cfg_s.ffr.Delta = DELTA;
cfg_s.ffr.adaptive = true;  cfg_s.ffr.adapt_dir = +1;
for i = 1:numel(aLo)
    for j = 1:numel(aHi)
        if aHi(j) < aLo(i), continue; end
        cfg_s.ffr.alpha_min = aLo(i);
        cfg_s.ffr.alpha_max = aHi(j);
        As = ffr_allocate(cfg_s, ctx, false);
        Fs = compute_sinr_ffr(cfg_s, ctx, As);
        Ks = compute_kpis(cfg_s, Fs, As);
        g1 = frontier_gain(sw.R_agg/1e3, sw.R_p5_edge, Ks.all.R_agg_Mbps/1e3, Ks.edge.R_p5_Mbps);
        g2 = frontier_gain(sw.R_p5,      sw.R_p5_edge, Ks.all.R_p5_Mbps,      Ks.edge.R_p5_Mbps);
        scan.gain_agg(i,j) = g1.gain_y;
        scan.gain_p5(i,j)  = g2.gain_y;
    end
end
[bestAgg, kA] = max(scan.gain_agg(:));  [iA,jA] = ind2sub(size(scan.gain_agg), kA);
[bestP5,  kP] = max(scan.gain_p5(:));   [iP,jP] = ind2sub(size(scan.gain_p5),  kP);
fprintf('  MEJOR ganancia sobre la frontera (R_agg, R_p5_borde): %+.3f Mbps  (anclas %.2f-%.2f)\n', ...
    bestAgg, aLo(iA), aHi(jA));
fprintf('  MEJOR ganancia sobre la frontera (R_p5,  R_p5_borde): %+.3f Mbps  (anclas %.2f-%.2f)\n', ...
    bestP5,  aLo(iP), aHi(jP));
scan.bestAgg = bestAgg;  scan.bestP5 = bestP5;
scan.bestAgg_anchors = [aLo(iA) aHi(jA)];
scan.bestP5_anchors  = [aLo(iP) aHi(jP)];
if bestAgg <= 0 && bestP5 <= 0
    fprintf(['  ==> NINGUN miembro de la familia supera la frontera estatica: el resultado\n' ...
             '      negativo es ESTRUCTURAL del escenario, no una mala eleccion de anclas.\n']);
end

%% 10. VEREDICTO H2 (honesto: se reporta tambien si NO mejora)
fprintf('\n============================== VEREDICTO H2 ==============================\n');
fprintf('Frontera estatica: %d valores de alpha fijo (%.2f a %.2f), Delta=%d, ventana completa.\n', ...
    nA, alphas(1), alphas(end), DELTA);

fprintf('\n(1) Plano pedido (R_agregado, R_p5_borde):\n');
if FR.degenPlane
    fprintf(['    PLANO DEGENERADO: %d punto no dominado de %d. R_agregado y R_p5_borde son\n' ...
             '    AMBOS decrecientes en alpha, luego alpha->0 domina y no hay frontera que\n' ...
             '    superar. La prueba es VACUA en este escenario (resultado en si mismo:\n' ...
             '    la sub-banda interior no compra agregado a ninguna elevacion).\n'], ...
        FR.nPareto, nA);
elseif FR.above
    fprintf(['    H2 SOPORTADA: a igual R_agregado (%.2f Gbps) la adaptativa da %+.2f Mbps\n' ...
             '    de R_p5 de BORDE (%.1f%% sobre el mejor alpha fijo).\n'], ...
        FR.x_q, FR.gain_y, 100*FR.gain_y_rel);
else
    fprintf('    NO mejora: deficit de %.2f Mbps de R_p5 de borde a igual agregado.\n', -FR.gain_y);
end

fprintf('\n(2) Plano de compromiso REAL (R_p5 global, R_p5 borde), %d puntos no dominados:\n', ...
    numel(GR.fx));
if GR.material
    fprintf(['    H2 SOPORTADA: a igual R_p5 GLOBAL (%.2f Mbps) la adaptativa da %+.3f Mbps\n' ...
             '    de R_p5 de BORDE (%+.1f%%) respecto al MEJOR alpha fijo equivalente,\n' ...
             '    por encima del umbral de materialidad (%.3f Mbps).\n'], ...
        GR.x_q, GR.gain_y, 100*GR.gain_y_rel, matThr);
elseif GR.above
    fprintf(['    H2 NO SOPORTADA: el punto adaptativo cae SOBRE la frontera estatica. La\n' ...
             '    ganancia (%+.3f Mbps, %+.1f%%) es POSITIVA pero NO MATERIAL: queda por\n' ...
             '    debajo del umbral %.3f Mbps (2%% de la frontera / medio paso del barrido),\n' ...
             '    es decir es indistinguible de la resolucion del propio barrido de alpha.\n' ...
             '    La adaptacion reproduce el compromiso de la mejor FFR estatica, no lo mejora.\n'], ...
        GR.gain_y, 100*GR.gain_y_rel, matThr);
else
    fprintf(['    H2 NO SOPORTADA: el punto adaptativo cae SOBRE o POR DEBAJO de la frontera\n' ...
             '    estatica (%+.3f Mbps de R_p5 de borde a igual R_p5 global). Con la regla\n' ...
             '    alpha-elevacion, la adaptacion NO compra un compromiso mejor que el mejor\n' ...
             '    alpha fijo. Es un resultado valido y se reporta como tal.\n'], GR.gain_y);
end

fprintf('\n(3) Por que (criterio analitico): la sub-banda INTERIOR renta solo si\n');
fprintf('    Delta*SE_centro/SE_borde > 1. Medido: %.3f (25-35 deg), %.3f (35-60), %.3f (60-90).\n', ...
    ratio_bin(1), ratio_bin(2), ratio_bin(3));
if all(ratio_bin < 1)
    fprintf(['    < 1 en TODOS los tramos: con %d haces co-canal la sub-banda interior\n' ...
             '    (reuso-1) nunca compensa su ancho, ni siquiera a elevacion alta. La premisa\n' ...
             '    de H2 ("en geometria favorable el interior rinde") NO se cumple en la\n' ...
             '    ventana de elevaciones simulada (max %.1f deg).\n'], ...
        BL.nBeams, max(ctx.elev_deg(ctx.cov)));
end

fprintf('\n(4) Contrastes de VALIDEZ (la regla no esta amanada):\n');
fprintf('    - sentido invertido (adapt_dir=-1): %+.3f Mbps sobre la frontera => la regla\n', GRi.gain_y);
fprintf('      con adapt_dir=+1 es %s que la invertida (corr(alpha,elev) = %+.3f).\n', ...
    ternary(GR.gain_y > GRi.gain_y,'MEJOR','PEOR'), corr_nan(alp_t, elev_t));
fprintf('    - caso degenerado alpha_min=alpha_max: %s (adaptativa == estatica, error 0).\n', ...
    ternary(degenOK,'PASA','FALLA'));
fprintf('    - familia de anclas: mejor ganancia alcanzable %+.3f Mbps (umbral material %.3f)\n', ...
    max(bestAgg,bestP5), matThr);
fprintf('      => %s.\n', ternary(max(bestAgg,bestP5) <= matThr, ...
    ['NINGUN miembro de la familia (alpha_min,alpha_max) supera la frontera de forma ' ...
     'material: el resultado negativo NO se debe a la eleccion de anclas 0.3/0.7'], ...
    'otras anclas darian margen material: conviene revisarlas'));

fprintf('\n(5) Lo que SI aporta la adaptativa (frente a la estatica alpha*=%.2f, no frente\n', alpha_star);
fprintf('    a la frontera): R_p5 de BORDE %.1f -> %.1f Mbps (x%.2f) al MISMO R_agregado\n', ...
    KK{iST}.edge.R_p5_Mbps, KK{iAD}.edge.R_p5_Mbps, ...
    KK{iAD}.edge.R_p5_Mbps/max(KK{iST}.edge.R_p5_Mbps,eps));
fprintf('    (%.2f vs %.2f Gbps), a costa de R_p5 global (%.1f -> %.1f Mbps).\n', ...
    KK{iAD}.all.R_agg_Mbps/1e3, KK{iST}.all.R_agg_Mbps/1e3, ...
    KK{iST}.all.R_p5_Mbps, KK{iAD}.all.R_p5_Mbps);
fprintf('==========================================================================\n');

%% 11. Figuras
%  Exportadas a PNG con save_fig (300 dpi, fuente unica) en FIGDIR, con nombres
%  DESCRIPTIVOS: son las figuras de uno de los dos experimentos cabecera.

% (b) LA GRAFICA CLAVE: frontera de compromiso + punto adaptativo
fA = figure('Name','E5 | Frontera de compromiso','Color','w');
plot(sw.R_agg/1e3, sw.R_p5_edge, 'o-', 'Color',[.45 .45 .45], 'LineWidth',1.4, ...
     'MarkerFaceColor','w', 'DisplayName','FFR estatica (barrido de \alpha fijo)'); hold on;
plot(FR.fx, FR.fy, 'k-', 'LineWidth',2.2, 'DisplayName','frontera de Pareto estatica');
for ia = 1:2:nA
    text(sw.R_agg(ia)/1e3, sw.R_p5_edge(ia), sprintf('  %.2f', alphas(ia)), ...
        'FontSize',7, 'Color',[.4 .4 .4]);
end
plot(KK{iR1}.all.R_agg_Mbps/1e3, KK{iR1}.edge.R_p5_Mbps, 'ks', ...
     'MarkerSize',10, 'MarkerFaceColor',[.8 .8 .8], 'DisplayName','reuso-1');
plot(KK{iST}.all.R_agg_Mbps/1e3, KK{iST}.edge.R_p5_Mbps, 'bd', ...
     'MarkerSize',10, 'MarkerFaceColor','b', 'DisplayName',sprintf('FFR estatica (\\alpha*=%.2f)',alpha_star));
plot(KK{iAD}.all.R_agg_Mbps/1e3, KK{iAD}.edge.R_p5_Mbps, 'rp', ...
     'MarkerSize',16, 'MarkerFaceColor','r', 'DisplayName','FFR ADAPTATIVA (H2)');
plot(KK{iIN}.all.R_agg_Mbps/1e3, KK{iIN}.edge.R_p5_Mbps, 'v', ...
     'Color',[.85 .5 .1], 'MarkerSize',9, 'MarkerFaceColor',[.85 .5 .1], ...
     'DisplayName','adaptativa con regla INVERTIDA');
% Segmento vertical que mide la ganancia sobre la frontera
if isfinite(FR.y_frontier)
    plot([FR.x_q FR.x_q], [FR.y_frontier KK{iAD}.edge.R_p5_Mbps], 'r--', ...
        'LineWidth',1.2, 'DisplayName',sprintf('ganancia = %+.2f Mbps', FR.gain_y));
end
grid on; xlabel('R_{agregado} [Gbps]'); ylabel('R_{p5} de los usuarios de BORDE [Mbps]');
title(sprintf('E5 | Frontera de compromiso de la FFR (\\Delta=%d): la adaptativa vs todo \\alpha fijo', DELTA));
legend('Location','best');
save_fig(fA, FIGDIR, 'e5_a_frontera_pareto_agregado_vs_borde');

% (b2) LA FRONTERA QUE SI ES UN COMPROMISO EN ESTE SISTEMA
fB = figure('Name','E5 | Frontera real (R_p5 global vs borde)','Color','w');
plot(sw.R_p5, sw.R_p5_edge, 'o-', 'Color',[.45 .45 .45], 'LineWidth',1.4, ...
     'MarkerFaceColor','w', 'DisplayName','FFR estatica (barrido de \alpha fijo)'); hold on;
plot(GR.fx, GR.fy, 'k-', 'LineWidth',2.2, 'DisplayName','frontera de Pareto estatica');
for ia = 1:2:nA
    text(sw.R_p5(ia), sw.R_p5_edge(ia), sprintf('  %.2f', alphas(ia)), ...
        'FontSize',7, 'Color',[.4 .4 .4]);
end
plot(KK{iR1}.all.R_p5_Mbps, KK{iR1}.edge.R_p5_Mbps, 'ks', ...
     'MarkerSize',10, 'MarkerFaceColor',[.8 .8 .8], 'DisplayName','reuso-1');
plot(KK{iST}.all.R_p5_Mbps, KK{iST}.edge.R_p5_Mbps, 'bd', ...
     'MarkerSize',10, 'MarkerFaceColor','b', 'DisplayName',sprintf('FFR estatica (\\alpha*=%.2f)',alpha_star));
plot(KK{iAD}.all.R_p5_Mbps, KK{iAD}.edge.R_p5_Mbps, 'rp', ...
     'MarkerSize',16, 'MarkerFaceColor','r', 'DisplayName','FFR ADAPTATIVA (H2)');
plot(KK{iIN}.all.R_p5_Mbps, KK{iIN}.edge.R_p5_Mbps, 'v', ...
     'Color',[.85 .5 .1], 'MarkerSize',9, 'MarkerFaceColor',[.85 .5 .1], ...
     'DisplayName','adaptativa con regla INVERTIDA');
if isfinite(GR.y_frontier)
    plot([GR.x_q GR.x_q], [GR.y_frontier KK{iAD}.edge.R_p5_Mbps], 'r--', ...
        'LineWidth',1.2, 'DisplayName',sprintf('ganancia = %+.3f Mbps', GR.gain_y));
end
grid on; xlabel('R_{p5} GLOBAL [Mbps]'); ylabel('R_{p5} de los usuarios de BORDE [Mbps]');
title(sprintf(['E5 | Compromiso REAL de la FFR (\\Delta=%d): proteccion de borde vs ' ...
               'percentil global'], DELTA));
legend('Location','best');
save_fig(fB, FIGDIR, 'e5_b_frontera_pareto_p5global_vs_borde');

% (c) alpha(t) y elevacion: la regla en accion
fC = figure('Name','E5 | alpha(t) y elevacion','Color','w');
th = (tvec - tvec(1))/60;
yyaxis left
plot(th, elev_t, 'LineWidth',1.6); ylabel('elevacion media del servidor [deg]');
yyaxis right
plot(th, alp_t, 'LineWidth',1.8); hold on;
plot(th, alpha_star*ones(size(th)), '--', 'LineWidth',1.4);
ylabel('\alpha (fraccion de banda del INTERIOR)'); ylim([0 1]);
grid on; xlabel('tiempo [min]');
title('E5 | La regla en accion: \alpha(t) sube con la elevacion (baja elev. -> banda al BORDE)');
legend({'elevacion','\alpha adaptativa', sprintf('\\alpha estatica = %.2f',alpha_star)}, ...
       'Location','best');
save_fig(fC, FIGDIR, 'e5_c_regla_alpha_vs_elevacion_temporal');

% (d) Desglose por tramo de elevacion (R_p5 de borde)
fD = figure('Name','E5 | Desglose por elevacion','Color','w');
bar(bin.R_p5_edge([iR1 iST iAD],:).');
set(gca,'XTickLabel',binName);
grid on; ylabel('R_{p5} de BORDE [Mbps]'); xlabel('tramo de elevacion del servidor');
title('E5 | Donde aporta la adaptacion: R_{p5} de borde por tramo de elevacion');
legend({cases{iR1}.name, cases{iST}.name, cases{iAD}.name}, 'Location','northwest');
save_fig(fD, FIGDIR, 'e5_d_rp5_borde_por_tramo_elevacion');

% (a) CDF de throughput de BORDE
fE = figure('Name','E5 | CDF throughput BORDE','Color','w');
for c = [iR1 iST iAD]
    cdf = KK{c}.edge.cdf_R;
    if all(isnan(cdf(:))), continue; end
    semilogx(cdf(:,1), cdf(:,2), 'LineWidth',1.7, 'DisplayName',cases{c}.name); hold on;
end
grid on; xlabel('Throughput R [Mbps]'); ylabel('CDF');
title('E5 | CDF de throughput de los usuarios de BORDE (p5 = metrica reina)');
legend('Location','southeast');
save_fig(fE, FIGDIR, 'e5_e_cdf_throughput_borde');

fprintf('\n[figuras] 5 PNG a 300 dpi en %s%s\n', FIGDIR, filesep);

%% 12. Guardar
ctxLite = rmfield(ctx, {'Grel','Glin'});

timing = struct('warmup_p618_s', tWarm, 'experiment_s', toc(tExp), 'total_s', toc(tRun));
fprintf('[tiempos] experimento %.1f s (+ %.1f s de warmup P.618) | total %.1f s\n', ...
    timing.experiment_s, timing.warmup_p618_s, timing.total_s);

save('e5_adaptive_results.mat', 'cfg','BL','ctxLite','cases','KK','KV','sw','alphas', ...
     'DELTA','alpha_star','alpha_starV','maskV','elMin_valid','FR','FRv','FRi', ...
     'GR','GRv','GRi','scan','ratio_bin','ratio_all', ...
     'bin','binName','elev_t','alp_t','degenOK','timing','MEM','BINFO');
fprintf('\nResultados guardados en e5_adaptive_results.mat\n');

%% ========================= funciones locales =========================
function printTable(cases, KK)
fprintf('%-16s %8s %9s %9s %11s %9s %7s %7s\n', ...
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
    fprintf('%-16s %8s %9.2f %9.1f %11.1f %9.2f %7.2f %7.3f\n', ...
        cases{c}.name, as, K.all.SINR_p5, K.all.R_p5_Mbps, K.edge.R_p5_Mbps, ...
        K.all.R_agg_Mbps/1e3, K.all.Pcov, K.all.jain);
end
fprintf('%s\n', repmat('-',1,80));
end

% -------------------------------------------------------------------------
function FR = frontier_gain(x, y, xq, yq)
%FRONTIER_GAIN  Posicion de un punto (xq,yq) respecto a la FRONTERA DE PARETO de
%   la nube (x,y), con el criterio "mas es mejor" en AMBOS ejes
%   (x = R_agregado, y = R_p5 de borde).
%
%   La frontera es el subconjunto de puntos NO DOMINADOS: (x_i,y_i) esta en la
%   frontera si ningun otro punto tiene a la vez x >= x_i e y >= y_i (con alguna
%   desigualdad estricta). Se interpola linealmente entre ellos.
%
%   Salidas:
%     .fx, .fy        frontera de Pareto (ordenada por x)
%     .y_frontier     y de la frontera evaluada en x = xq
%     .x_frontier     x de la frontera evaluada en y = yq
%     .gain_y         yq - y_frontier   (>0 => el punto esta POR ENCIMA)
%     .gain_x         xq - x_frontier
%     .above          logico: el punto mejora estrictamente la frontera
x = x(:).';  y = y(:).';
ok = isfinite(x) & isfinite(y);
x = x(ok);   y = y(ok);

% Puntos no dominados
keep = true(1, numel(x));
for i = 1:numel(x)
    dom = (x >= x(i)) & (y >= y(i)) & ((x > x(i)) | (y > y(i)));
    if any(dom), keep(i) = false; end
end
fx = x(keep);  fy = y(keep);
[fx, is] = sort(fx);  fy = fy(is);
% Duplicados en x (misma abscisa): quedarse con el mejor y
[fxu, ~, gidx] = unique(fx);
fyu = accumarray(gidx(:), fy(:), [], @max).';
FR.fx = fxu;  FR.fy = fyu;

FR.x_q = xq;  FR.y_q = yq;
FR.y_frontier = interp_clamped(fxu, fyu, xq);
FR.x_frontier = interp_clamped(fliplr(fyu), fliplr(fxu), yq);   % fy es decreciente en x
FR.gain_y = yq - FR.y_frontier;
FR.gain_x = xq - FR.x_frontier;
FR.gain_y_rel = FR.gain_y / max(abs(FR.y_frontier), eps);
FR.gain_x_rel = FR.gain_x / max(abs(FR.x_frontier), eps);
% "Por encima" en sentido estricto: el punto no esta dominado por NINGUN alpha fijo
% y ademas mejora la interpolacion de la frontera.
FR.dominated = any((x >= xq) & (y >= yq) & ((x > xq) | (y > yq)));
FR.above     = ~FR.dominated && (FR.gain_y > 0);
end

% -------------------------------------------------------------------------
function yq = interp_clamped(xx, yy, xq)
%INTERP_CLAMPED  Interpolacion lineal con extrapolacion PLANA en los extremos
%   (fuera del rango cubierto por la frontera no se extrapola: se satura).
if numel(xx) < 2
    if isempty(xx), yq = NaN; else, yq = yy(1); end
    return;
end
if xq <= xx(1),   yq = yy(1);   return; end
if xq >= xx(end), yq = yy(end); return; end
yq = interp1(xx, yy, xq, 'linear');
end

% -------------------------------------------------------------------------
function report_frontier(label, FR)
fprintf('  %-32s  agregado %.2f Gbps | R_p5_borde %.1f Mbps | frontera %.1f | %s %+0.2f Mbps\n', ...
    label, FR.x_q, FR.y_q, FR.y_frontier, ...
    ternary(FR.above,'ENCIMA','debajo'), FR.gain_y);
end

% -------------------------------------------------------------------------
function report_frontier2(label, FR)
fprintf('  %-30s %8.2f %10.2f %14.2f %+10.3f  %s\n', ...
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
