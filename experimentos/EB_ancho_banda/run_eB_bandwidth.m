%% RUN_EB_BANDWIDTH  EB: sensibilidad de la FFR al ANCHO DE BANDA DEL SISTEMA.
%
%  PREGUNTA DEL TUTOR
%  ------------------
%  "¿Influye el ancho de banda B en la ventaja de la FFR? Con mas BW caben mas
%   usuarios y la FFR es menos critica; con menos BW caben menos y repartir bien
%   importa mas." Exigencia adicional: el NUMERO DE USUARIOS debe mantenerse
%   constante entre puntos (aqui se verifica y se aborta si no).
%
%  PREDICCION A FALSAR (escrita ANTES de medir; NO se fuerza en el codigo)
%  ----------------------------------------------------------------------
%  En compute_sinr_ffr los cuatro terminos de la SINR llevan el MISMO factor de
%  ancho asignado B_user (que es proporcional a B por la particion de
%  ffr_allocate, B_sub = alpha*B o (1-alpha)*B/Delta, repartida entre n0*f_clase
%  usuarios):
%
%      C       = Cpsd  + 10log10(B_user/1e6) + Grel_propio       (PSD constante)
%      I_intra = Cpsd  + 10log10(B_user/1e6) + 10log10(sum G co-canal)
%      I_inter = Ipsd  + 10log10(B_user/1e6) + 10log10(f)
%      N       = -228.6 + 10log10(Tsys) + 10log10(B_user)
%
%  y ademas Ipsd_inter es INVARIANTE a B por construccion (ffr_context hace
%  Ipsd = INTF.I_dBW - 10log10(B_MHz), y compute_interference habia formado
%  I_dBW con EIRP = EIRPdens + 10log10(B_MHz): el factor se cancela). Por tanto
%
%      SINR = C/(N + I_intra + I_inter)   es INVARIANTE a B,
%      R    = B_user * log2(1+SINR)       es EXACTAMENTE LINEAL en B.
%
%  Se espera pues un RESULTADO NULO en SINR, cobertura, veredictos y
%  crit_ratio_ffr, y escalado exactamente proporcional en R. EL OBJETIVO DE ESTE
%  EXPERIMENTO ES MEDIRLO, NO ASUMIRLO: si alguna comprobacion sale distinta de
%  cero mas alla del ruido de coma flotante, es un BUG y hay que investigarlo.
%
%  Nota de coma flotante: los cinco B del barrido son 250/4, 250/2, 250, 250*2 y
%  250*4, es decir escalados por potencias de 2 (EXACTOS en IEEE-754), luego los
%  anchos B_user son exactamente proporcionales. Lo unico que puede introducir
%  ruido son las evaluaciones de log10/10^ del camino dB, del orden de 1e-13 dB.
%  Por eso las desviaciones se REPORTAN CON SU VALOR y no se redondean a "0".
%
%  ARQUITECTURA (politica del proyecto: motor compartido + runner delgado)
%  ----------------------------------------------------------------------
%  NO hay fisica en este fichero. Cada punto del barrido lo ejecuta
%  run_one_density (funcion pura de punto: pipeline completo, copia local de cfg,
%  assert_cfg_coherent, los 6 esquemas sobre el MISMO ctx) y el bucle es
%  run_sweep_points. B se cambia UNICAMENTE por cfg.radio.B_MHz a traves de
%  pt.cfg_over (merge recursivo con merge_cfg); radio_band_Hz es la fuente unica
%  de la conversion a Hz y NO existe ningun cfg.radio.B_Hz que pueda quedar
%  petrificado (auditoria, hallazgo G3 -- que es EXACTAMENTE la trampa que este
%  barrido habria pisado).
%
%  IMPORTANTE (honestidad del experimento): se ejecuta el PIPELINE COMPLETO en
%  cada punto -- geometria, radioenlace, interferencia y ffr_context incluidos --
%  aunque analiticamente ninguno de ellos dependa de B. Reutilizar un unico ctx y
%  reevaluar solo ffr_allocate/compute_sinr_ffr seria 5 veces mas rapido, pero
%  daria por supuesta la mitad de la invariancia que el experimento debe MEDIR.
%
%  Escenario: el de cabecera (Ku real, Walker 53:1584/72/1, 91 haces) con el
%  perfil de COSTE convergido de la FASE A (step_km=4 -> M=317, Nt=61), que es
%  ademas el de E3b: eso da un ANCLA DE REGRESION exacta (ver seccion 5).
%
%  Uso: ejecutar en la raiz del proyecto. Resultados en eB_bandwidth.mat,
%  figuras en figs_eB/. No modifica config_default ni ningun .mat publicado.

clear; clc; close all;

tRun   = tic;
FIGDIR = 'figs_eB';

fprintf('=======================================================================\n');
fprintf('  EB: SENSIBILIDAD AL ANCHO DE BANDA DEL SISTEMA\n');
fprintf('=======================================================================\n\n');

%% ------------------------------------------------------------------------
%  1. CONFIGURACION BASE -- todo CONGELADO salvo cfg.radio.B_MHz
%  ------------------------------------------------------------------------
%  Se parte de config_default() y se sobrescribe una COPIA, igual que hacen
%  config_calib y config_oneweb. config_default.m NO se toca.
cfg = config_default();

% Constelacion REAL de cabecera (la de config_default; se explicita por claridad)
cfg.constellations(1).T = 1584;
cfg.constellations(1).P = 72;

% Perfil de COSTE convergido de la FASE A. La rejilla NO depende de B, luego el
% numero de usuarios M es el mismo en los cinco puntos por construccion; aun asi
% se VERIFICA a posteriori y se aborta si no (exigencia del tutor).
cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 4;          % FASE A -> M = 317 usuarios
cfg.time.dt          = 120;        % FASE A -> Nt = 61 instantes
cfg.time.duration    = 7200;       % 2 h

% Umbral centro/borde por CUANTIL al 50% (mismo criterio que E2/E5/E3b). Es lo
% que mantiene la poblacion de borde IDENTICA entre esquemas y, aqui, entre
% puntos de B: el cuantil se toma sobre SINR_ref, que es invariante a B.
cfg.ffr.tau_mode = 'quantile';
cfg.ffr.tau_q    = 50;

% Cluster de 91 haces (guarda 2.07 anillos con el HPBW real de 2.08 deg)
cfg.beams.nRings = 5;
cfg.beams.nBeams = [];

% Troceado temporal. Se pide 18 -- el mismo que uso de hecho la tanda de E3b que
% sirve de ancla de regresion -- para que cualquier residuo de la regresion sea
% atribuible SOLO a B. El troceado es EXACTO de todas formas
% (test_timeblock_invariance y test_ctx_blocked: max|dif| = 0), asi que si el
% estimador de memoria obliga a bajarlo, el resultado no cambia.
timeBlockReq = 18;                                                  % <-- DECISION
cfg.compute.timeBlock = timeBlockReq;

%% ------------------------------------------------------------------------
%  2. EJE DEL BARRIDO Y ESQUEMAS
%  ------------------------------------------------------------------------
%  B_MHz = 250 es el ancho de canal real de la solicitud FCC de SpaceX y es el
%  punto de referencia (todo el corpus del TFG esta calculado ahi). Los otros
%  cuatro son x1/4, x1/2, x2 y x4 de ese valor: cubren un rango de x16 y son
%  potencias de 2 exactas del punto de referencia, lo que hace que el escalado
%  de anchos sea exacto en coma flotante y el residuo medible sea solo el de
%  log10/10^.
Blist  = [62.5 125 250 500 1000];                                   % MHz
nBpt   = numel(Blist);
Bref   = 250;
iRef   = find(Blist == Bref, 1);
if isempty(iRef)
    error('run_eB_bandwidth:noRef', 'El barrido debe incluir el punto de referencia B = %g MHz.', Bref);
end

alphaStar = 0.5;                   % mismo alpha* que E3b/E3 (barrido de E2: 0.5-0.6)

cases = { ...
    struct('scheme','reuse1','Delta',1,'alpha',NaN,      'adaptive',false,'name','reuse1'), ...
    struct('scheme','reuseD','Delta',3,'alpha',NaN,      'adaptive',false,'name','reuseD(3)'), ...
    struct('scheme','reuseD','Delta',4,'alpha',NaN,      'adaptive',false,'name','reuseD(4)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',alphaStar,'adaptive',false,'name','ffr(3,a*)'), ...
    struct('scheme','ffr',   'Delta',4,'alpha',alphaStar,'adaptive',false,'name','ffr(4,a*)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',alphaStar,'adaptive',true, 'name','ffr-adapt(3)') };
nC = numel(cases);
nm = cellfun(@(c) c.name, cases, 'UniformOutput', false);

% Indices utiles para las comparaciones de cabecera
iR1 = 1;  iF3 = 4;  iF4 = 5;  iA3 = 6;

fprintf('Perfil: T=%d P=%d | rejilla %.0f km paso %.0f km | dt=%.0f s | HPBW=%.2f deg | nRings=%d\n', ...
    cfg.constellations(1).T, cfg.constellations(1).P, cfg.ground.radius_km, ...
    cfg.ground.step_km, cfg.time.dt, cfg.radio.beamwidth3dB_deg, cfg.beams.nRings);
fprintf('Eje EB: B = [%s] MHz (referencia %g) | alpha* = %.2f | tau cuantil %g%% | sched=%s\n', ...
    num2str(Blist), Bref, alphaStar, cfg.ffr.tau_q, cfg.ffr.sched);
fprintf('Esquemas (%d): %s\n\n', nC, strjoin(nm, ' | '));

%% ------------------------------------------------------------------------
%  3. PUNTOS DEL BARRIDO Y PLAN DE MEMORIA
%  ------------------------------------------------------------------------
%  B se cambia SOLO por cfg.radio.B_MHz. merge_cfg es recursivo, luego el resto
%  de cfg.radio (freq_GHz, EIRPdensity, Grx_max, GT/Tsys...) queda intacto; y
%  run_one_density llama a assert_cfg_coherent justo despues de aplicar el
%  override, antes de simular nada.
pts = cell(1, nBpt);
for i = 1:nBpt
    pts{i} = struct( ...
        'cfg_over', struct('radio', struct('B_MHz', Blist(i))), ...
        'cases',    {cases}, ...
        'label',    sprintf('B%g', Blist(i)));
end

% B_MHz NO forma parte de la clave de la cache P.618
% (freq|lat|lon|minElev|p618_avail), luego los 5 puntos COMPARTEN una sola tabla.
fprintf('--------- PLAN DE MEMORIA ---------\n');
MEM = estimate_sweep_memory(cfg, pts, 1, struct('verbose',true));
if ~isempty(MEM.timeBlock_reco) && ~isnan(MEM.timeBlock_reco)
    fprintf(['\n[AVISO] El estimador reduce cfg.compute.timeBlock de %d a %d por RAM libre.\n' ...
             '        El troceado es EXACTO (test_timeblock_invariance / test_ctx_blocked:\n' ...
             '        max|dif| = 0), asi que la regresion contra el ancla sigue debiendo dar 0\n' ...
             '        aunque el ancla se calculara con otro tamano de bloque.\n'], ...
        cfg.compute.timeBlock, MEM.timeBlock_reco);
    cfg.compute.timeBlock = MEM.timeBlock_reco;
    MEM = estimate_sweep_memory(cfg, pts, 1, struct('verbose',false));
end
timeBlockUsed = cfg.compute.timeBlock;

% SECUENCIAL a proposito: el limite de esta maquina es la RAM (cada worker
% mantiene su propia geometria M*N*Ntb), y con una sola tabla P.618 el paralelo
% no compra nada aqui.
nWorkers = 0;                                                       % <-- DECISION

%% ------------------------------------------------------------------------
%  4. BARRIDO
%  ------------------------------------------------------------------------
%  Warmup EXPLICITO de la tabla P.618 fuera del cronometro (coste fijo de la
%  configuracion atmosferica, ~80 s, cacheado en un persistent; test_atm_cache
%  verifica que la version cacheada es identica bit a bit).
fprintf('\n[warmup] construyendo la tabla ITU-R P.618 fuera del cronometro...\n');
tW = tic;  atm_loss_dB(cfg, 45);  tWarm = toc(tW);
fprintf('[warmup] %.1f s (coste fijo, se descuenta)\n', tWarm);

tExp = tic;
[res, INFO] = run_sweep_points(cfg, pts, nWorkers, struct('verbose',true));
tSweep = toc(tExp);

%% 4b. INVARIANTES DEL ESCENARIO (exigencia del tutor: misma poblacion)
%  Si alguno de estos cambiara, el barrido no estaria aislando B y cualquier
%  diferencia posterior seria incomparable. Se ABORTA, no se avisa.
M0 = res{iRef}.M;  Nt0 = res{iRef}.Nt;  N0 = res{iRef}.N;  nB0 = res{iRef}.nBeams;
for i = 1:nBpt
    if res{i}.M ~= M0 || res{i}.Nt ~= Nt0 || res{i}.N ~= N0 || res{i}.nBeams ~= nB0
        error('run_eB_bandwidth:scenarioDrift', ...
            ['El escenario cambia con B (B=%g MHz: M=%d Nt=%d N=%d nBeams=%d frente a ' ...
             'M=%d Nt=%d N=%d nBeams=%d en la referencia): el barrido NO esta aislando el ' ...
             'ancho de banda.'], Blist(i), res{i}.M, res{i}.Nt, res{i}.N, res{i}.nBeams, ...
             M0, Nt0, N0, nB0);
    end
end
fprintf(['\n[check] escenario IDENTICO en los %d puntos: M=%d usuarios, Nt=%d instantes, ' ...
         'N=%d satelites, %d haces.\n'], nBpt, M0, Nt0, N0, nB0);

%% ------------------------------------------------------------------------
%  5. EXTRACCION DE METRICAS  [nBpt x nC]
%  ------------------------------------------------------------------------
z = @() nan(nBpt, nC);
Q = struct( ...
    'SINR_edge_p5', z(), 'SINR_p5_glob', z(), 'Pcov_edge',   z(), ...
    'R_p5_edge',    z(), 'R_p5_glob',    z(), 'R_agg_Gbps',  z(), ...
    'SE_mean',      z(), 'crit_ratio',   z(), 'n0',          z(), ...
    'B_user_MHz',   z(), 'penalty_mean', z(), 'SINR_c_p5',   z(), ...
    'centerFrac',   z(), 'n_edge',       z(), 'alpha_mean',  z(), ...
    'SE_edge_p5',   z(), 'covFrac',      z());
verdict = cell(nBpt, nC);

for i = 1:nBpt
    O = res{i};
    for c = 1:nC
        K = O.K{c};
        Q.SINR_edge_p5(i,c) = K.viab.SINR_edge_p5;
        Q.SINR_p5_glob(i,c) = K.all.SINR_p5;
        Q.Pcov_edge(i,c)    = K.edge.Pcov;
        Q.R_p5_edge(i,c)    = K.edge.R_p5_Mbps;
        Q.R_p5_glob(i,c)    = K.all.R_p5_Mbps;
        Q.R_agg_Gbps(i,c)   = K.all.R_agg_Mbps / 1e3;
        Q.SE_mean(i,c)      = K.all.SE_mean_bpsHz;
        Q.crit_ratio(i,c)   = O.diag.crit_ratio_ffr(c);
        Q.n0(i,c)           = O.cases{c}.load_n0_mean;
        Q.B_user_MHz(i,c)   = K.B_user_MHz_mean;
        Q.penalty_mean(i,c) = K.penalty_mean;
        Q.SINR_c_p5(i,c)    = K.center.SINR_p5;
        Q.centerFrac(i,c)   = K.centerFrac;
        Q.n_edge(i,c)       = K.viab.n_edge;
        Q.alpha_mean(i,c)   = O.cases{c}.alpha_mean;
        Q.SE_edge_p5(i,c)   = K.viab.SE_edge_p5;
        Q.covFrac(i,c)      = K.covFrac;
        verdict{i,c}        = K.viab.verdict;
    end
end

%% ------------------------------------------------------------------------
%  6. TABLA POR PUNTO DE B
%  ------------------------------------------------------------------------
for i = 1:nBpt
    fprintf('\n#################  B = %g MHz  (x%.4g respecto de %g MHz)  #################\n', ...
        Blist(i), Blist(i)/Bref, Bref);
    fprintf('%-14s %9s %9s %8s %9s %9s %9s %8s %9s %7s %9s  %s\n', ...
        'esquema','SINRe_p5','SINRp5_g','Pcob_B','R_p5_bor','R_p5_glo','R_agreg','SE_med', ...
        'critFFR','n0','B_usr','VEREDICTO');
    fprintf('%-14s %9s %9s %8s %9s %9s %9s %8s %9s %7s %9s\n', ...
        '','[dB]','[dB]','','[Mbps]','[Mbps]','[Gbps]','[b/s/Hz]','','[u/celda]','[MHz]');
    fprintf('%s\n', repmat('-',1,120));
    for c = 1:nC
        if isnan(Q.crit_ratio(i,c)), crTxt = '      -'; else, crTxt = sprintf('%7.4f', Q.crit_ratio(i,c)); end
        fprintf('%-14s %9.4f %9.4f %8.4f %9.3f %9.3f %9.4f %8.4f %9s %7.3f %9.4f  %s\n', ...
            nm{c}, Q.SINR_edge_p5(i,c), Q.SINR_p5_glob(i,c), Q.Pcov_edge(i,c), ...
            Q.R_p5_edge(i,c), Q.R_p5_glob(i,c), Q.R_agg_Gbps(i,c), Q.SE_mean(i,c), ...
            crTxt, Q.n0(i,c), Q.B_user_MHz(i,c), upper(verdict{i,c}));
    end
    fprintf('%s\n', repmat('-',1,120));
    fprintf('penalizacion media (C/N - SINR) [dB]: ');
    penPairs = [nm; num2cell(Q.penalty_mean(i,:))];
    fprintf('%s=%.3f  ', penPairs{:});
    fprintf('\n');
end

%% ------------------------------------------------------------------------
%  7. COMPROBACIONES (a) - (d)
%  ------------------------------------------------------------------------
fprintf('\n\n=======================================================================\n');
fprintf('  COMPROBACIONES DEL RESULTADO NULO (referencia: B = %g MHz)\n', Bref);
fprintf('=======================================================================\n');

CHK = struct();

% --- (a) SINR de borde ---------------------------------------------------
dS = abs(Q.SINR_edge_p5 - Q.SINR_edge_p5(iRef,:));
CHK.max_dSINR_edge_p5 = max(dS(:));
[~, lin] = max(dS(:));  [ia, ca] = ind2sub(size(dS), lin);
fprintf('\n(a) SINR_edge_p5: max|SINR(B) - SINR(%g)| = %.6e dB', Bref, CHK.max_dSINR_edge_p5);
fprintf('  (peor: %s a B=%g MHz)\n', nm{ca}, Blist(ia));
fprintf('    SINR_p5 GLOBAL: max|dif| = %.6e dB\n', max(abs(Q.SINR_p5_glob - Q.SINR_p5_glob(iRef,:)),[],'all'));
CHK.max_dSINR_p5_glob = max(abs(Q.SINR_p5_glob - Q.SINR_p5_glob(iRef,:)),[],'all');

% --- (b) cobertura, criterio de H2 y veredictos --------------------------
CHK.max_dPcov_edge  = max(abs(Q.Pcov_edge  - Q.Pcov_edge(iRef,:)),  [], 'all');
crd = abs(Q.crit_ratio - Q.crit_ratio(iRef,:));
CHK.max_dCritRatio  = max(crd(~isnan(crd)));
if isempty(CHK.max_dCritRatio), CHK.max_dCritRatio = NaN; end
CHK.nVerdictChanges = 0;
vlist = {};
for i = 1:nBpt
    for c = 1:nC
        if ~strcmp(verdict{i,c}, verdict{iRef,c})
            CHK.nVerdictChanges = CHK.nVerdictChanges + 1;
            vlist{end+1} = sprintf('%s@B=%g: %s (ref %s)', nm{c}, Blist(i), ...
                verdict{i,c}, verdict{iRef,c}); %#ok<SAGROW>
        end
    end
end
fprintf('\n(b) Pcov de BORDE:      max|dif| = %.6e\n', CHK.max_dPcov_edge);
fprintf('    crit_ratio_ffr:     max|dif| = %.6e   (solo esquemas ffr; NaN en el resto)\n', CHK.max_dCritRatio);
fprintf('    VEREDICTOS:         %d cambios de %d (%d puntos x %d esquemas)\n', ...
    CHK.nVerdictChanges, nBpt*nC, nBpt, nC);
if ~isempty(vlist), fprintf('      %s\n', vlist{:}); end
% Poblacion de borde: si el cuantil tau se moviera con B, la clasificacion
% cambiaria y las comparaciones dejarian de ser sobre la misma poblacion.
CHK.max_dNedge      = max(abs(Q.n_edge     - Q.n_edge(iRef,:)),     [], 'all');
CHK.max_dCenterFrac = max(abs(Q.centerFrac - Q.centerFrac(iRef,:)), [], 'all');
fprintf('    poblacion de BORDE: max|dif| n_edge = %g muestras | max|dif| centerFrac = %.6e\n', ...
    CHK.max_dNedge, CHK.max_dCenterFrac);

% --- (c) linealidad del throughput ---------------------------------------
ratioB   = Blist(:) ./ Bref;                          % [nBpt x 1] escalado nominal
ratioObs = Q.R_p5_edge ./ Q.R_p5_edge(iRef,:);        % observado
devAbs   = abs(ratioObs - ratioB);
devRel   = devAbs ./ ratioB;
CHK.max_dev_ratio_abs = max(devAbs(:));
CHK.max_dev_ratio_rel = max(devRel(:));
ratioObsG   = Q.R_p5_glob ./ Q.R_p5_glob(iRef,:);
CHK.max_dev_ratio_glob_rel = max(abs(ratioObsG - ratioB) ./ ratioB, [], 'all');
ratioObsA   = Q.R_agg_Gbps ./ Q.R_agg_Gbps(iRef,:);
CHK.max_dev_ratio_agg_rel  = max(abs(ratioObsA - ratioB) ./ ratioB, [], 'all');

fprintf('\n(c) LINEALIDAD de R con B (cociente observado frente a B/%g):\n', Bref);
fprintf('%10s %10s', 'B [MHz]', 'B/Bref');
for c = 1:nC, fprintf('%14s', nm{c}); end
fprintf('\n%s\n', repmat('-',1,20+14*nC));
for i = 1:nBpt
    fprintf('%10g %10.4f', Blist(i), ratioB(i));
    for c = 1:nC, fprintf('%14.9f', ratioObs(i,c)); end
    fprintf('\n');
end
fprintf('%s\n', repmat('-',1,20+14*nC));
fprintf('    R_p5 de BORDE : max|cociente - B/Bref| = %.6e  (relativo %.6e)\n', ...
    CHK.max_dev_ratio_abs, CHK.max_dev_ratio_rel);
fprintf('    R_p5 GLOBAL   : desviacion relativa max = %.6e\n', CHK.max_dev_ratio_glob_rel);
fprintf('    R_agregado    : desviacion relativa max = %.6e\n', CHK.max_dev_ratio_agg_rel);
% El ancho por usuario debe escalar EXACTAMENTE (potencias de 2)
buRatio = Q.B_user_MHz ./ Q.B_user_MHz(iRef,:);
CHK.max_dev_Buser_rel = max(abs(buRatio - ratioB) ./ ratioB, [], 'all');
fprintf('    B_user medio  : desviacion relativa max = %.6e (escalado exacto esperado)\n', ...
    CHK.max_dev_Buser_rel);

% --- (d) GANANCIA DE LA FFR SOBRE REUSO-1 (la metrica que responde al tutor) ---
gain_dB = Q.SINR_edge_p5 - Q.SINR_edge_p5(:,iR1);      % [nBpt x nC]
gain_x  = Q.R_p5_edge    ./ Q.R_p5_edge(:,iR1);
CHK.gain_dB = gain_dB;
CHK.gain_x  = gain_x;
CHK.max_spread_gain_dB = max(max(gain_dB,[],1) - min(gain_dB,[],1));
CHK.max_spread_gain_x  = max(max(gain_x, [],1) - min(gain_x, [],1));

fprintf('\n(d) GANANCIA DE LA FFR SOBRE REUSO-1 por punto de B:\n');
fprintf('%10s', 'B [MHz]');
for c = 2:nC, fprintf('%13s', [nm{c} ' dB']); end
for c = 2:nC, fprintf('%13s', [nm{c} ' x']); end
fprintf('\n%s\n', repmat('-',1,10+26*(nC-1)));
for i = 1:nBpt
    fprintf('%10g', Blist(i));
    for c = 2:nC, fprintf('%13.4f', gain_dB(i,c)); end
    for c = 2:nC, fprintf('%13.4f', gain_x(i,c));  end
    fprintf('\n');
end
fprintf('%s\n', repmat('-',1,10+26*(nC-1)));
fprintf('    RECORRIDO de la ganancia en todo el eje de B (x16):\n');
for c = 2:nC
    fprintf('      %-14s  dB: %.4f -> %.4f (recorrido %.6e dB) | x: %.4f -> %.4f (recorrido %.6e)\n', ...
        nm{c}, min(gain_dB(:,c)), max(gain_dB(:,c)), max(gain_dB(:,c))-min(gain_dB(:,c)), ...
        min(gain_x(:,c)),  max(gain_x(:,c)),  max(gain_x(:,c))-min(gain_x(:,c)));
end

%% ------------------------------------------------------------------------
%  8. REGRESION contra el resultado PUBLICADO (ancla exacta)
%  ------------------------------------------------------------------------
%  Esta es la comprobacion que demuestra que el resultado nulo NO procede de que
%  el barrido no este cambiando nada: si el punto B=250 reproduce bit a bit un
%  resultado ya publicado, el motor esta haciendo lo de siempre, y las otras
%  cuatro columnas se han obtenido con el MISMO codigo cambiando solo B_MHz.
%
%  ANCLA: minelev_results.mat (E3b), indice R{1}{5} = densidad 'denso' (T=1584),
%  minElev = 25 deg. Ese punto usa EXACTAMENTE esta configuracion -- perfil FASE A
%  (step_km=4, dt=120), T=1584/72, nRings=5, tau cuantil 50%, alpha*=0.5 y los
%  MISMOS 6 esquemas -- y su minElev de 25 deg coincide con el de config_default,
%  luego el override de E3b es un no-op sobre este perfil.
fprintf('\n\n=======================================================================\n');
fprintf('  REGRESION DEL PUNTO B = %g MHz CONTRA RESULTADO PUBLICADO\n', Bref);
fprintf('=======================================================================\n');

REG = struct('done',false,'anchor','minelev_results.mat R{1}{5} (E3b denso, minElev=25)', ...
             'max_abs_diff',NaN,'nFieldsCompared',0,'notes',{{}},'exact',false);
if ~exist('minelev_results.mat','file')
    fprintf(['\n[REGRESION NO REALIZADA] No existe minelev_results.mat en la carpeta. ' ...
             'NO se declara superada.\n']);
else
    Sa = load('minelev_results.mat','R','minElevList','densList','cases','alphaStar','cfg');
    ok = true;  why = {};
    if Sa.alphaStar ~= alphaStar,                     ok=false; why{end+1}='alpha* distinto'; end
    if Sa.cfg.ground.step_km ~= cfg.ground.step_km,   ok=false; why{end+1}='step_km distinto'; end
    if Sa.cfg.time.dt ~= cfg.time.dt,                 ok=false; why{end+1}='dt distinto'; end
    if Sa.cfg.beams.nRings ~= cfg.beams.nRings,       ok=false; why{end+1}='nRings distinto'; end
    if Sa.cfg.radio.B_MHz ~= Bref,                    ok=false; why{end+1}='B_MHz del ancla distinto'; end
    if numel(Sa.cases) ~= nC,                         ok=false; why{end+1}='numero de esquemas distinto'; end
    iAnc = find(Sa.minElevList == cfg.geom.minElev, 1);
    if isempty(iAnc),                                 ok=false; why{end+1}='el ancla no tiene minElev=25'; end

    if ~ok
        fprintf(['\n[REGRESION NO REALIZADA] La configuracion del ancla NO coincide (%s). ' ...
                 'NO se declara superada.\n'], strjoin(why,', '));
        REG.notes = why;
    else
        REF = Sa.R{1}{iAnc};
        fprintf('\nAncla: %s | label=%s | T=%d M=%d Nt=%d nBeams=%d | timeBlock ancla=%d, aqui=%d\n', ...
            REG.anchor, REF.label, REF.T, REF.M, REF.Nt, REF.nBeams, ...
            Sa.cfg.compute.timeBlock, timeBlockUsed);
        CMP = struct('max',0,'n',0,'notes',{{}});
        for c = 1:nC
            if ~strcmp(REF.cases{c}.name, nm{c})
                CMP.notes{end+1} = sprintf('esquema %d: nombre "%s" vs "%s"', c, REF.cases{c}.name, nm{c});
            end
            CMP = cmp_rec(res{iRef}.K{c}, REF.K{c}, sprintf('K{%d}', c), CMP);
        end
        % diag: se compara la interseccion; los campos NUEVOS se listan aparte.
        CMP = cmp_rec(res{iRef}.diag, REF.diag, 'diag', CMP);
        fNew = setdiff(fieldnames(res{iRef}.cases{1}), fieldnames(REF.cases{1}));

        REG.done            = true;
        REG.max_abs_diff    = CMP.max;
        REG.nFieldsCompared = CMP.n;
        REG.notes           = CMP.notes;
        REG.exact           = (CMP.max == 0) && isempty(CMP.notes);

        fprintf('Campos numericos comparados: %d | max|dif| = %.6e\n', CMP.n, CMP.max);
        if ~isempty(fNew)
            fprintf(['Campos NUEVOS en cases{} respecto del ancla (anadidos en esta tanda, ' ...
                     'informativos, no entran en ningun calculo): %s\n'], strjoin(fNew', ', '));
        end
        if ~isempty(CMP.notes)
            fprintf('DISCREPANCIAS ESTRUCTURALES:\n');  fprintf('  %s\n', CMP.notes{:});
        end
        if REG.exact
            fprintf('==> REGRESION SUPERADA: max|dif| = 0 EXACTO contra el resultado publicado.\n');
        else
            fprintf(['==> REGRESION **NO** EXACTA (max|dif| = %.6e). NO se declara superada: ' ...
                     'hay que investigar el origen antes de dar por bueno el resultado nulo.\n'], CMP.max);
        end
    end
end

%% ------------------------------------------------------------------------
%  9. PARTE 2 -- LA HIPOTESIS DEL TUTOR, HECHA METRICA: CAPACIDAD EN USUARIOS
%  ------------------------------------------------------------------------
%  El simulador NO tiene control de admision ni modelo de demanda (limitacion
%  declarada), luego "cuantos usuarios caben" NO es una salida del motor. Se
%  DERIVA como post-proceso sobre KPIs ya calculados, sin tocar nada.
%
%  DERIVACION (auditable; es analitica, no simulada)
%  -------------------------------------------------
%  Con la convencion de planificacion sched='share' de ffr_allocate, un usuario
%  de la clase X en una celda con n usuarios recibe
%
%       B_user(X) = B_sub(X) / (n * f_X)
%
%  donde f_X es la fraccion de usuarios de esa clase en la celda y B_sub(X) el
%  ancho de la sub-banda que usa la clase. Exigir la tasa objetivo R_min al
%  usuario del PERCENTIL 5 de la clase (la metrica reina del TFG) da
%
%       B_sub(X)/(n*f_X) * SE_p5(X) >= R_min   =>   n <= B_sub(X)*SE_p5(X)/(R_min*f_X)
%
%  con SE_p5(X) = log2(1 + 10^(SINR_p5(X)/10)) [bit/s/Hz]. La celda debe cumplirlo
%  para TODAS sus clases a la vez, luego
%
%       n_max = min_X  B_sub(X)*SE_p5(X) / (R_min*f_X)
%
%  Instanciacion por esquema (f_c = fraccion de centro medida, f_e = 1 - f_c):
%    reuse1 : todas las clases comparten la banda PLENA        -> B_sub = B,       f_X = 1
%    reuseD : todas las clases comparten una sub-banda de color-> B_sub = B/Delta, f_X = 1
%    ffr    : centro -> B_sub = alpha*B,            f = f_c
%             borde  -> B_sub = (1-alpha)*B/Delta,  f = f_e
%  (en reuse1/reuseD la clase que manda es siempre la de BORDE, porque su SE_p5 es
%   menor y comparte el mismo ancho; en ffr puede mandar cualquiera de las dos, y
%   eso es justo el compromiso que la FFR ajusta con alpha.)
%
%  SUPUESTOS DECLARADOS: (i) la SINR no depende de la carga (cierto en este
%  modelo: n0 solo divide el ancho, no entra en la SINR); (ii) se usa el
%  percentil 5 de cada clase, no la media, luego n_max es la carga que sostiene
%  el objetivo para el 95% de los usuarios de cada clase; (iii) en el esquema
%  ADAPTATIVO alpha varia por instante y aqui se usa su media temporal, asi que
%  su n_max es una aproximacion (se marca en la tabla).
fprintf('\n\n=======================================================================\n');
fprintf('  PARTE 2: CAPACIDAD EN USUARIOS POR CELDA (post-proceso analitico)\n');
fprintf('=======================================================================\n');

RminList = [1 5 10 25] * 1e6;                 % bit/s objetivo por usuario   % <-- DECISION
nR       = numel(RminList);
nmax     = nan(nBpt, nC, nR);
bindCls  = cell(nBpt, nC, nR);                % clase que limita ('centro'/'borde')

for i = 1:nBpt
    B_Hz = Blist(i) * 1e6;
    for c = 1:nC
        fC = Q.centerFrac(i,c);  fE = 1 - fC;
        SEc = log2(1 + 10^(Q.SINR_c_p5(i,c)/10));
        SEe = log2(1 + 10^(Q.SINR_edge_p5(i,c)/10));
        D   = res{i}.cases{c}.Delta;
        switch lower(res{i}.cases{c}.scheme)
            case 'reuse1', Bc_ = B_Hz;                  Be_ = B_Hz;                  fc_ = 1;  fe_ = 1;
            case 'reused', Bc_ = B_Hz/D;                Be_ = B_Hz/D;                fc_ = 1;  fe_ = 1;
            case 'ffr'
                a   = Q.alpha_mean(i,c);
                Bc_ = a*B_Hz;   Be_ = (1-a)*B_Hz/D;     fc_ = fC; fe_ = fE;
        end
        for r = 1:nR
            nc = Bc_*SEc / (RminList(r)*max(fc_,eps));
            ne = Be_*SEe / (RminList(r)*max(fe_,eps));
            [nmax(i,c,r), w] = min([nc ne]);
            wlbl = {'centro','borde'};  bindCls{i,c,r} = wlbl{w};
        end
    end
end

ratio_nmax = nmax ./ nmax(:,iR1,:);           % frente a reuso-1, por punto de B

for r = 1:nR
    fprintf('\n--- R_min = %.0f Mbps por usuario ---\n', RminList(r)/1e6);
    fprintf('%10s', 'B [MHz]');
    for c = 1:nC, fprintf('%16s', nm{c}); end
    fprintf('\n%10s', '');
    for c = 1:nC, fprintf('%16s', 'n_max (clase)'); end
    fprintf('\n%s\n', repmat('-',1,10+16*nC));
    for i = 1:nBpt
        fprintf('%10g', Blist(i));
        for c = 1:nC
            fprintf('%11.2f(%s)', nmax(i,c,r), bindCls{i,c,r}(1));
        end
        fprintf('\n');
    end
    fprintf('%s\n', repmat('-',1,10+16*nC));
    fprintf('%10s', 'COCIENTE');
    for c = 1:nC, fprintf('%16.6f', ratio_nmax(iRef,c,r)); end
    fprintf('   <- n_max(esquema)/n_max(reuso-1) en B=%g\n', Bref);
    sp = max(ratio_nmax(:,:,r),[],1) - min(ratio_nmax(:,:,r),[],1);
    fprintf('%10s', 'recorrido');
    for c = 1:nC, fprintf('%16.3e', sp(c)); end
    fprintf('   <- variacion del cociente en TODO el eje de B (x16)\n');
end

CHK.max_spread_ratio_nmax = max(abs(ratio_nmax - ratio_nmax(iRef,:,:)), [], 'all');
fprintf(['\n[PARTE 2] max|n_max(esquema)/n_max(reuso-1) en B  -  el mismo cociente en B=%g| ' ...
         '= %.6e\n'], Bref, CHK.max_spread_ratio_nmax);
fprintf(['          Si ese recorrido es ~0, la ventaja de la FFR medida EN USUARIOS ADMITIDOS\n' ...
         '          es la MISMA con 62.5 MHz que con 1000 MHz, y la hipotesis "con mas BW la\n' ...
         '          FFR es menos critica" queda refutada DENTRO de este modelo.\n']);

%% ------------------------------------------------------------------------
%  10. PARTE 3 -- DONDE B SI IMPORTARIA Y NO ESTA MODELADO (DIAGNOSTICO)
%  ------------------------------------------------------------------------
%  ESTO NO ES UN RESULTADO DE SIMULACION. El motor trata la banda como un
%  continuo (ffr_allocate parte B en numeros reales), luego ningun B pequeno
%  puede penalizar a la particion de la FFR por granularidad. En un sistema 5G NR
%  real la banda se asigna en BLOQUES DE RECURSO (PRB) y ESA cuantizacion es el
%  unico mecanismo identificado por el que un B pequeno perjudicaria de verdad a
%  la FFR: al partir B en Delta sub-bandas, cada una debe seguir siendo un numero
%  util de PRB.
%
%  NUMEROLOGIA DECLARADA, y su salvedad (cuestion abierta del trabajo):
%  la banda Ku de 10.7-12.7 GHz NO es una banda NR estandarizada -- queda por
%  encima de FR1 (<=7.125 GHz) y por debajo de FR2 (>=24.25 GHz), y las bandas
%  NTN de la Rel-17 son de L/S. Por tanto la numerologia AQUI ES UNA ELECCION
%  DECLARADA, no una lectura de la norma: se tabulan dos casos, SCS = 30 kHz
%  (estilo FR1) y SCS = 120 kHz (estilo FR2), con
%       1 PRB = 12 subportadoras  ->  360 kHz (30 kHz) / 1.44 MHz (120 kHz)
%  y el limite de 275 PRB por portadora NR (TS 38.211). El umbral de "numero
%  razonable de PRB" se declara abajo.
fprintf('\n\n=======================================================================\n');
fprintf('  PARTE 3: DIAGNOSTICO DE CUANTIZACION EN PRB (NO es simulacion)\n');
fprintf('=======================================================================\n');

SCSlist   = [30 120] * 1e3;        % Hz   numerologias contrastadas          % <-- DECISION
NPRB_max  = 275;                   % PRB  maximo por portadora NR (TS 38.211)
PRB_min   = 11;                    % PRB  umbral declarado de "asignacion razonable"  % <-- DECISION
% PROVENIENCIA de PRB_min: el menor ancho de canal de FR1 (5 MHz) corresponde a
% 11 PRB con SCS = 30 kHz en la tabla de configuracion de transmision de
% TS 38.104 (Tabla 5.3.2-1). PENDIENTE: verificar la cifra en la norma antes de
% citarla en la memoria; el diagnostico se puede releer con cualquier otro umbral
% porque la tabla imprime los PRB en bruto.
DeltaList = [3 4];

fprintf(['\nUmbral declarado: %d PRB (= 5 MHz con SCS=30 kHz, el menor canal de FR1).\n' ...
         'Carga medida n0 = %.3f usuarios/celda, fraccion de borde f_e = %.3f (tau cuantil %g%%).\n' ...
         'alpha = %.2f.  "!" marca sub-banda o asignacion por usuario por debajo del umbral.\n'], ...
    PRB_min, Q.n0(iRef,iF3), 1-Q.centerFrac(iRef,iF3), cfg.ffr.tau_q, alphaStar);

PRB = struct('SCS_Hz',SCSlist,'Delta',DeltaList,'PRB_min',PRB_min,'NPRB_max',NPRB_max, ...
             'is_diagnostic',true, ...
             'note','Diagnostico analitico: el motor trata la banda como CONTINUA. No es salida de la simulacion.');
PRB.Bc_MHz      = nan(nBpt, numel(DeltaList));
PRB.Be_MHz      = nan(nBpt, numel(DeltaList));
PRB.Buser_e_MHz = nan(nBpt, numel(DeltaList));
PRB.nprb_sub_e  = nan(nBpt, numel(DeltaList), numel(SCSlist));
PRB.nprb_usr_e  = nan(nBpt, numel(DeltaList), numel(SCSlist));
PRB.flag        = false(nBpt, numel(DeltaList), numel(SCSlist));

for s = 1:numel(SCSlist)
    prbHz = 12 * SCSlist(s);
    fprintf('\n--- SCS = %g kHz  (1 PRB = %.3f MHz) ---\n', SCSlist(s)/1e3, prbHz/1e6);
    fprintf('%9s %7s %10s %10s %12s %11s %11s %s\n', ...
        'B [MHz]','Delta','Bc [MHz]','Be [MHz]','B_usr_e[MHz]','PRB sub-e','PRB usr-e','aviso');
    fprintf('%s\n', repmat('-',1,92));
    for i = 1:nBpt
        for d = 1:numel(DeltaList)
            D  = DeltaList(d);
            Bc = alphaStar*Blist(i);
            Be = (1-alphaStar)*Blist(i)/D;
            cIdx = find(cellfun(@(x) strcmpi(x.scheme,'ffr') && x.Delta==D && ~x.adaptive, res{i}.cases), 1);
            fE   = 1 - Q.centerFrac(i,cIdx);
            n0i  = Q.n0(i,cIdx);
            Bue  = Be / (n0i*fE);
            np_s = floor(Be *1e6/prbHz);
            np_u = floor(Bue*1e6/prbHz);
            fl   = (np_u < PRB_min) || (np_s < PRB_min);
            msg  = '';
            if np_u < PRB_min, msg = [msg sprintf('! usuario %d<%d PRB ', np_u, PRB_min)]; end
            if np_s < PRB_min, msg = [msg sprintf('! sub-banda %d<%d PRB ', np_s, PRB_min)]; end
            if np_s > NPRB_max, msg = [msg sprintf('(sub-banda %d PRB > %d: exige agregacion) ', np_s, NPRB_max)]; end
            fprintf('%9g %7d %10.4f %10.4f %12.4f %11d %11d %s\n', ...
                Blist(i), D, Bc, Be, Bue, np_s, np_u, msg);
            PRB.Bc_MHz(i,d) = Bc;  PRB.Be_MHz(i,d) = Be;  PRB.Buser_e_MHz(i,d) = Bue;
            PRB.nprb_sub_e(i,d,s) = np_s;  PRB.nprb_usr_e(i,d,s) = np_u;
            PRB.flag(i,d,s) = fl;
        end
    end
    fprintf('%s\n', repmat('-',1,92));
end
fprintf(['\n[ETIQUETA] Las cifras de esta seccion son un DIAGNOSTICO ANALITICO, no un\n' ...
         '           resultado de la simulacion: el modelo trata la banda como continua y\n' ...
         '           por tanto NO reproduce ninguna penalizacion por granularidad. La\n' ...
         '           cuantizacion en PRB es el unico mecanismo identificado por el que un B\n' ...
         '           pequeno castigaria de verdad a la particion de la FFR, y queda como\n' ...
         '           TRABAJO FUTURO (junto con las tablas MODCOD).\n']);

%% ------------------------------------------------------------------------
%  11. FIGURAS
%  ------------------------------------------------------------------------
mk = {'o','s','d','^','v','p'};
ls = {'-','-','-','--','--',':'};
% Ancho generoso: con la leyenda 'eastoutside' el eje se estrecha y un titulo
% largo se sale del lienzo -- save_fig recorta al contenido y lo cortaria.
FIGPOS = [80 80 1000 520];

% (1) EL RESULTADO NULO: SINR_edge_p5 vs B, plana
f1 = figure('Name','SINR de borde vs ancho de banda','Color','w','Position',FIGPOS);
for c = 1:nC
    semilogx(Blist, Q.SINR_edge_p5(:,c), [ls{c} mk{c}], 'LineWidth',1.6, ...
        'MarkerSize',7, 'DisplayName',nm{c}); hold on;
end
yline(cfg.viab.gamma_th_dB,    'k--','LineWidth',1.4, ...
    'DisplayName',sprintf('\\gamma_{th} = %.1f dB (umbral)', cfg.viab.gamma_th_dB));
yline(cfg.viab.gamma_floor_dB, 'r:', 'LineWidth',1.4, ...
    'DisplayName',sprintf('\\gamma_{floor} = %.1f dB (QPSK 5G NR)', cfg.viab.gamma_floor_dB));
grid on; xlabel('Ancho de banda del sistema B [MHz]'); ylabel('SINR_{edge,p5} [dB]');
xticks(Blist); xticklabels(string(Blist)); xlim([min(Blist)*0.8 max(Blist)*1.25]);
title({'Resultado nulo: SINR de borde invariante a B', ...
       sprintf('x%g de B mueve el KPI %.2e dB', max(Blist)/min(Blist), CHK.max_dSINR_edge_p5)});
legend('Location','eastoutside');
save_fig(f1, FIGDIR, 'eB_a_sinr_borde_vs_ancho_de_banda');

% (2) LINEALIDAD: R_p5 de borde vs B en log-log con recta de pendiente 1
f2 = figure('Name','Throughput de borde vs ancho de banda','Color','w','Position',FIGPOS);
for c = 1:nC
    loglog(Blist, Q.R_p5_edge(:,c), [ls{c} mk{c}], 'LineWidth',1.6, ...
        'MarkerSize',7, 'DisplayName',nm{c}); hold on;
end
refY = Q.R_p5_edge(iRef,iR1) * (Blist/Bref);
loglog(Blist, refY, 'k-', 'LineWidth',1.0, 'DisplayName','pendiente 1 (R \propto B)');
grid on; xlabel('Ancho de banda del sistema B [MHz]'); ylabel('R_{p5} de BORDE [Mbps]');
xticks(Blist); xticklabels(string(Blist));
title({'R_{p5} de borde: exactamente proporcional a B', ...
       sprintf('desviacion relativa max. frente a pendiente 1: %.2e', CHK.max_dev_ratio_rel)});
legend('Location','eastoutside');
save_fig(f2, FIGDIR, 'eB_b_throughput_borde_vs_ancho_loglog');

% (3) LA RESPUESTA AL TUTOR: ganancia de la FFR sobre reuso-1 vs B
f3 = figure('Name','Ganancia de la FFR vs ancho de banda','Color','w','Position',FIGPOS);
yyaxis left
plot(log10(Blist), gain_dB(:,iF3), '-o','LineWidth',1.8,'MarkerSize',7, ...
    'DisplayName','\Delta=3: ganancia [dB]'); hold on;
plot(log10(Blist), gain_dB(:,iF4), '-s','LineWidth',1.8,'MarkerSize',7, ...
    'DisplayName','\Delta=4: ganancia [dB]');
ylabel('Ganancia de SINR_{edge,p5} sobre reuso-1 [dB]');
yyaxis right
plot(log10(Blist), gain_x(:,iF3), '--^','LineWidth',1.8,'MarkerSize',7, ...
    'DisplayName','\Delta=3: R_{p5} borde [x]');
plot(log10(Blist), gain_x(:,iF4), '--v','LineWidth',1.8,'MarkerSize',7, ...
    'DisplayName','\Delta=4: R_{p5} borde [x]');
ylabel('R_{p5} de borde frente a reuso-1 [veces]');
grid on; xlabel('Ancho de banda del sistema B [MHz]');
xticks(log10(Blist)); xticklabels(string(Blist));
title({'La ventaja de la FFR no depende de B', ...
       sprintf('recorrido en x%g de B: %.2e dB y %.2e veces', ...
       max(Blist)/min(Blist), CHK.max_spread_gain_dB, CHK.max_spread_gain_x)});
legend('Location','eastoutside');
save_fig(f3, FIGDIR, 'eB_c_ganancia_ffr_vs_ancho_de_banda');

% (4) PARTE 2: n_max(FFR)/n_max(reuso-1) vs B, una curva por R_min
f4 = figure('Name','Usuarios admitidos: FFR frente a reuso-1','Color','w','Position',FIGPOS);
for r = 1:nR
    semilogx(Blist, squeeze(ratio_nmax(:,iF3,r)), '-o','LineWidth',1.6,'MarkerSize',7, ...
        'DisplayName',sprintf('\\Delta=3, R_{min} = %g Mbps', RminList(r)/1e6)); hold on;
end
for r = 1:nR
    semilogx(Blist, squeeze(ratio_nmax(:,iF4,r)), '--s','LineWidth',1.3,'MarkerSize',6, ...
        'DisplayName',sprintf('\\Delta=4, R_{min} = %g Mbps', RminList(r)/1e6));
end
yline(1,'k:','LineWidth',1.4,'DisplayName','paridad con reuso-1');
grid on; xlabel('Ancho de banda del sistema B [MHz]');
ylabel('n_{max}(FFR) / n_{max}(reuso-1)');
xticks(Blist); xticklabels(string(Blist)); xlim([min(Blist)*0.8 max(Blist)*1.25]);
title({'Usuarios admitidos: la ventaja de la FFR es plana en B', ...
       sprintf('recorrido max. del cociente en x%g de B: %.2e', ...
       max(Blist)/min(Blist), CHK.max_spread_ratio_nmax)});
legend('Location','eastoutside');
save_fig(f4, FIGDIR, 'eB_d_usuarios_admitidos_ffr_vs_reuso1');

fprintf('\n[figuras] 4 PNG a 300 dpi en %s%s\n', FIGDIR, filesep);

%% ------------------------------------------------------------------------
%  12. VEREDICTO DEL EXPERIMENTO
%  ------------------------------------------------------------------------
tTot = toc(tRun);
timing = struct('warmup_p618_s',tWarm,'sweep_s',tSweep,'total_s',tTot);

% Umbral de MATERIALIDAD declarado: por debajo de esto la variacion no es
% distinguible del ruido de coma flotante del camino en dB (log10/10^ sobre
% magnitudes de ~1e2 dB dan ~1e-13; se toma 1e-9 como margen holgado).
TOL_dB = 1e-9;                                                      % <-- DECISION
hypothesisHolds = (CHK.max_spread_gain_dB > TOL_dB) || (CHK.max_spread_gain_x > 1e-6);

fprintf('\n\n=======================================================================\n');
fprintf('  VEREDICTO DEL EXPERIMENTO EB\n');
fprintf('=======================================================================\n');
fprintf('Eje: B = [%s] MHz (x%g) | %d esquemas | M=%d usuarios CONSTANTE | Nt=%d | T=%d\n', ...
    num2str(Blist), max(Blist)/min(Blist), nC, M0, Nt0, res{iRef}.T);
fprintf('\nDESVIACIONES MAXIMAS frente al punto de referencia B = %g MHz:\n', Bref);
fprintf('  (a) SINR_edge_p5 .................. %.6e dB\n', CHK.max_dSINR_edge_p5);
fprintf('      SINR_p5 global ................ %.6e dB\n', CHK.max_dSINR_p5_glob);
fprintf('  (b) Pcov de borde ................. %.6e\n',    CHK.max_dPcov_edge);
fprintf('      crit_ratio_ffr ................ %.6e\n',    CHK.max_dCritRatio);
fprintf('      veredictos que cambian ........ %d de %d\n', CHK.nVerdictChanges, nBpt*nC);
fprintf('      poblacion de borde (n_edge) ... %g muestras\n', CHK.max_dNedge);
fprintf('  (c) R_p5 borde frente a B/Bref .... %.6e absoluto (%.6e relativo)\n', ...
    CHK.max_dev_ratio_abs, CHK.max_dev_ratio_rel);
fprintf('      R_agregado (relativo) ......... %.6e\n', CHK.max_dev_ratio_agg_rel);
fprintf('  (d) ganancia FFR sobre reuso-1 .... recorrido %.6e dB / %.6e veces\n', ...
    CHK.max_spread_gain_dB, CHK.max_spread_gain_x);
fprintf('      usuarios admitidos (Parte 2) .. recorrido %.6e\n', CHK.max_spread_ratio_nmax);

fprintf('\nREGRESION (B = %g MHz contra %s):\n', Bref, REG.anchor);
if ~REG.done
    fprintf('  NO REALIZADA -> NO se declara superada. Motivo: %s\n', ...
        strjoin([{'ver seccion 8'} REG.notes], '; '));
elseif REG.exact
    fprintf('  SUPERADA: max|dif| = %.6e sobre %d campos numericos (0 EXACTO).\n', ...
        REG.max_abs_diff, REG.nFieldsCompared);
else
    fprintf('  NO SUPERADA: max|dif| = %.6e sobre %d campos numericos.\n', ...
        REG.max_abs_diff, REG.nFieldsCompared);
end

fprintf('\nHIPOTESIS DEL TUTOR ("con mas BW la FFR es menos critica"):\n');
if hypothesisHolds
    fprintf(['  SE SOSTIENE en este modelo: la ganancia de la FFR sobre reuso-1 varia %.4f dB\n' ...
             '  y %.4f veces al recorrer B en un factor %g. Investigar el mecanismo.\n'], ...
        CHK.max_spread_gain_dB, CHK.max_spread_gain_x, max(Blist)/min(Blist));
else
    fprintf(['  NO SE SOSTIENE en este modelo. Con la poblacion de usuarios CONSTANTE, la\n' ...
             '  ventaja de la FFR sobre reuso-1 es la MISMA con %g MHz que con %g MHz:\n' ...
             '    - SINR_edge_p5: %+.2f dB (Delta=3) y %+.2f dB (Delta=4), con un recorrido\n' ...
             '      de %.2e dB en todo el eje (x%g de ancho de banda);\n' ...
             '    - R_p5 de borde: x%.3f (Delta=3) y x%.3f (Delta=4), recorrido %.2e;\n' ...
             '    - usuarios admitidos por celda a igual R_min: cociente FFR/reuso-1\n' ...
             '      x%.3f (Delta=3), recorrido %.2e.\n' ...
             '  CAUSA: B se cancela en la SINR (entra por igual en C, I_intra, I_inter y N)\n' ...
             '  y solo escala la tasa; luego duplicar B duplica el caudal de TODOS los\n' ...
             '  esquemas por igual y NO cambia el reparto entre centro y borde, que es lo\n' ...
             '  unico que la FFR decide. "Caben mas usuarios" es cierto -- n_max crece\n' ...
             '  proporcionalmente a B -- pero crece IGUAL en reuso-1 y en FFR.\n'], ...
        min(Blist), max(Blist), gain_dB(iRef,iF3), gain_dB(iRef,iF4), ...
        CHK.max_spread_gain_dB, max(Blist)/min(Blist), ...
        gain_x(iRef,iF3), gain_x(iRef,iF4), CHK.max_spread_gain_x, ...
        ratio_nmax(iRef,iF3,1), CHK.max_spread_ratio_nmax);
    fprintf(['  DONDE SI IMPORTARIA (no modelado, ver Parte 3): la cuantizacion en PRB. El\n' ...
             '  motor parte la banda como un continuo; con B pequeno y Delta grande la\n' ...
             '  sub-banda de borde y, sobre todo, la asignacion POR USUARIO caen a pocos PRB.\n']);
end
fprintf('\n[tiempos] barrido %.1f s (+ %.1f s de tabla P.618) | total %.1f s\n', ...
    tSweep, tWarm, tTot);
fprintf('=======================================================================\n');

%% ------------------------------------------------------------------------
%  13. GUARDADO
%  ------------------------------------------------------------------------
save('eB_bandwidth.mat', 'res','INFO','cfg','cases','nm','Blist','Bref','iRef', ...
     'alphaStar','Q','verdict','CHK','REG','nmax','ratio_nmax','RminList','bindCls', ...
     'PRB','MEM','timing','timeBlockReq','timeBlockUsed','-v7.3');
fprintf('\nResultados en eB_bandwidth.mat | figuras en %s%s\n', FIGDIR, filesep);

%% ========================================================================
%  FUNCIONES LOCALES
%  ========================================================================
function CMP = cmp_rec(A, B, path, CMP)
%CMP_REC  Comparacion recursiva campo a campo de dos estructuras de resultados.
%   Acumula en CMP.max el maximo |dif| numerico, en CMP.n el numero de elementos
%   numericos comparados y en CMP.notes las discrepancias ESTRUCTURALES (campos
%   que faltan, tamanos distintos, cadenas distintas). Las discrepancias
%   estructurales NO se convierten en un numero: se listan, porque un campo que
%   falta no es "una diferencia de 0".
if isstruct(A) && isstruct(B)
    fa = fieldnames(A);  fb = fieldnames(B);
    only = setdiff(fa, fb);
    if ~isempty(only)
        CMP.notes{end+1} = sprintf('%s: campos solo en el nuevo: %s', path, strjoin(only', ','));
    end
    for i = 1:numel(fb)
        k = fb{i};
        if ~isfield(A, k)
            CMP.notes{end+1} = sprintf('%s.%s: falta en el resultado nuevo', path, k);
            continue;
        end
        CMP = cmp_rec(A.(k), B.(k), [path '.' k], CMP);
    end
elseif iscell(A) && iscell(B)
    if ~isequal(size(A), size(B))
        CMP.notes{end+1} = sprintf('%s: tamano de cell distinto (%s vs %s)', path, ...
            mat2str(size(A)), mat2str(size(B)));
        return;
    end
    for i = 1:numel(A)
        CMP = cmp_rec(A{i}, B{i}, sprintf('%s{%d}', path, i), CMP);
    end
elseif ischar(A) || ischar(B) || isstring(A) || isstring(B)
    if ~isequal(A, B)
        CMP.notes{end+1} = sprintf('%s: texto distinto ("%s" vs "%s")', path, string(A), string(B));
    end
elseif (isnumeric(A) || islogical(A)) && (isnumeric(B) || islogical(B))
    if ~isequal(size(A), size(B))
        CMP.notes{end+1} = sprintf('%s: tamano distinto (%s vs %s)', path, ...
            mat2str(size(A)), mat2str(size(B)));
        return;
    end
    a = double(A(:));  b = double(B(:));
    both = isnan(a) & isnan(b);                 % NaN == NaN cuenta como igual
    if any(isnan(a) ~= isnan(b))
        CMP.notes{end+1} = sprintf('%s: patron de NaN distinto', path);
    end
    d = abs(a(~both) - b(~both));
    if ~isempty(d), CMP.max = max(CMP.max, max(d)); end
    CMP.n = CMP.n + numel(a);
else
    if ~isequaln(A, B)
        CMP.notes{end+1} = sprintf('%s: tipos o valores no comparables', path);
    end
end
end
