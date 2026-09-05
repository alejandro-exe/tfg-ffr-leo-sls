%% RUN_E4_INTERCONSTELLATION  E4: interferencia INTER-CONSTELACION (eventos in-line).
%
%  QUE CONSTELACIONES ENTRAN: CRITERIO DE CO-CANAL EN BANDA Ku
%  ------------------------------------------------------------
%  Este experimento mide interferencia CO-CANAL sobre NUESTRO downlink Ku (12 GHz).
%  Por tanto solo puede figurar como interferente un operador cuyo ENLACE DE USUARIO
%  DE BAJADA comparta esa banda. Los dos que lo hacen son:
%
%     Starlink (propio) : Ku 10.7-12.7 GHz al terminal de usuario  [SpaceX FCC]
%     OneWeb            : Ku 10.7-12.7 GHz al terminal de usuario  [OneWeb FCC]
%
%  KUIPER QUEDA FUERA, y no es un recorte de alcance sino una CORRECCION FISICA:
%  el enlace de usuario de bajada de Kuiper es de banda **Ka** (17.7-20.2 GHz), luego
%  NO es co-canal con nuestro downlink Ku y no puede generar interferencia co-canal
%  sobre el. Una version anterior de este experimento lo incluia como tercer operador
%  "en la misma banda"; era incorrecto y todo resultado que dependiera de el (en
%  particular la elasticidad del numero de operadores) queda RETIRADO.
%
%  CONSECUENCIA: el eje "numero de operadores" tiene DOS puntos reales (propia sola /
%  propia + OneWeb). Con dos puntos no se puede ajustar una elasticidad, y no se
%  intenta: el hallazgo se sostiene sobre el CONTRASTE Starlink-OneWeb, que es dato
%  real (ver TAREA 4).
%
%  FISICA (por que este eje es DISTINTO de E3 y E3b)
%  --------------------------------------------------
%  E3 demostro que densificar no satura porque la reticula co-canal es INVARIANTE DE
%  ESCALA: con s = sqrt(3)*sin(HPBW/2) el vecino del 1er anillo cae siempre en
%  u = 1.61634*sqrt(3) = 2.7996 -> -10.67 dB, sea cual sea el ancho de haz. Pero esa
%  auto-similitud protege UNICAMENTE a los haces de la PROPIA reticula: los haces de
%  otro operador NO pertenecen a ella y pueden caer a cualquier angulo.
%
%  Aun asi la interferencia inter-constelacion MEDIA sera pequena, y por una razon ya
%  medida en este proyecto: el terminal VSAT (ITU-R S.1428, term_gain_dB) discrimina
%  35-48 dB contra cualquier satelite fuera de eje, sea de quien sea. Lo que importa
%  no es la media sino la COLA:
%
%      EVENTOS IN-LINE. Cuando un satelite de otro operador se acerca angularmente a
%      la linea de vista del servidor (psi -> 0), el terminal DEJA DE DISCRIMINAR
%      (term_gain_dB(0) = Grx_max) y la interferencia se dispara varias decenas de dB
%      de golpe. Es exactamente el fenomeno que regulan los limites EPFD de la UIT en
%      la coordinacion NGSO-NGSO.
%
%  DOS ESCALAS DE TIEMPO (decision metodologica CRITICA)
%  -----------------------------------------------------
%  Un LEO barre ~0.5-1 deg/s de separacion angular vista desde tierra, luego un
%  evento con psi < 2 deg dura SEGUNDOS. Con el dt = 120 s del pipeline de KPIs los
%  eventos se pierden casi por completo: se veria una probabilidad y una duracion sin
%  sentido. Por eso E4 trabaja en dos niveles:
%
%    NIVEL FINO (TAREA 2, el nucleo): punto de referencia UNICO (M = 1) y dt = 1 s
%      sobre la ventana completa. Sin rejilla de usuarios el coste es N*Nt y cabe de
%      sobra. Da la CDF de psi_min, P(psi_min < X), la DURACION de los eventos y la
%      SINR dentro/fuera, con las MISMAS formulas de compute_interference (psi por la
%      esferica, term_gain_dB S.1428, beam_gain_dB del haz ajeno, satPointing).
%
%    NIVEL KPI (TAREAS 3-4): el pipeline completo con la rejilla, dt = 120 s y
%      pt.constellations = 1 / 2 operadores. Da SINR_edge_p5, p1 y VERDICT.
%      **Este nivel INFRAMUESTREA los eventos in-line**, luego su impacto en los KPIs
%      es una COTA INFERIOR; la magnitud real de la cola la da el nivel fino. Se dice
%      explicitamente en las tablas para que no se lea al reves.
%
%  SALVEDAD DE MODELADO QUE HAY QUE DECLARAR
%  ------------------------------------------
%  compute_sinr_ffr reescala la interferencia INTER-satelite por f = 1/Delta en las
%  sub-bandas de borde, suponiendo que el interferente usa NUESTRO coloreado y solo
%  coincide con probabilidad 1/Delta. Para un operador AJENO SIN COORDINAR eso no se
%  cumple: ocupa toda la banda y el factor correcto seria f = 1. La aproximacion hace
%  que la FFR parezca algo mejor de lo que es frente a interferencia ajena. Es
%  inmaterial mientras I_inter << N (que es lo que ocurre de media), pero NO tiene por
%  que serlo durante un evento in-line: por eso el nivel fino reporta I_ajena/N
%  directamente, sin pasar por ese factor, y asi se puede juzgar el sesgo.
%
%  Uso: ejecutar en la raiz del proyecto. Resultados en interconstellation_results.mat,
%  figuras en figs_e4/. No modifica la fisica del motor.

clear; clc; close all;

fprintf('=======================================================================\n');
fprintf('  E4: INTERFERENCIA INTER-CONSTELACION (eventos in-line)\n');
fprintf('=======================================================================\n\n');

%% ------------------------------------------------------------------------
%  1. TAREA 1: ESCENARIO MULTI-OPERADOR
%  ------------------------------------------------------------------------
cfg = config_default();

cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 4;
cfg.time.dt          = 120;
cfg.time.duration    = 7200;
cfg.ffr.tau_mode     = 'quantile';
cfg.ffr.tau_q        = 50;
cfg.beams.nRings     = 5;
cfg.beams.nBeams     = [];
cfg.compute.timeBlock = 20;

% CRITICO en multi-operador: el usuario es cliente del operador 1 y solo puede
% servirse de EL. Sin esto associate_serving elegiria el mejor satelite de
% CUALQUIER constelacion -- el usuario se conectaria a la competencia -- y ademas
% ese servidor apuntaria celdas Earth-fixed proyectadas para h = 550 km desde
% 1200 km, hundiendo la SINR por un motivo que NO es interferencia. Verificado: sin
% esta linea el pipeline atribuia 5-9 dB de perdida a la inter-constelacion cuando
% el nivel fino mide I_ajena/N <= -13.6 dB (imposible).            % <-- DECISION
cfg.geom.serving_cid = 1;

% Los DOS operadores que comparten banda Ku de usuario (ver cabecera: Kuiper opera en
% Ka y por eso NO entra). El 1 es el PROPIO (el que sirve al usuario); el 2 es ajeno y
% emite en la MISMA banda -> peor caso sin coordinacion espectral.
%
% ONEWEB con datos REALES de su solicitud (antes: 648 sat / 12 planos, estimado):
%   720 satelites = 18 planos x 40 sat, 1200 km, 87.9 deg, constelacion polar 'star'
%   [OneWeb Schedule S Technical Report, p.2]
%   EIRPdensity = 10.58 dBW/MHz  (-13.4 dBW/4 kHz)  [Technical Narrative, p.20 Sec. A.12]
%   16 haces Ku, huella 1140x1140 km, Gtx 24.5 dBi, 250 MHz por canal
%
% MASCARA DE ELEVACION POR OPERADOR (campo minElev, ver build_constellation):
%   OneWeb NO opera por debajo de 55 deg (55 nominal global; 45-50 en bajas latitudes)
%   [Technical Narrative, p.13 Sec. A.4], frente a los 25 deg de Starlink. Sin esta
%   mascara entrarian como interferentes satelites de OneWeb a elevaciones a las que
%   NO transmiten, y P(evento in-line) saldria SOBREESTIMADA.
OPS = struct( ...
    'name',    {'Starlink-S1','OneWeb'}, ...
    'pattern', {'delta','star'}, ...
    'h',       {550, 1200}, ...
    'inc',     {53, 87.9}, ...
    'T',       {1584, 720}, ...
    'P',       {72, 18}, ...
    'F',       {1, 1}, ...
    'Om0',     {0, 0}, ...
    'M0',      {0, 0}, ...
    'minElev', {25, 55});

dt_fine = 1;                      % s, nivel fino (eventos in-line)      % <-- DECISION
psiTh   = [1 2 5 10];             % deg, umbrales de evento in-line
nOpsList = [1 2];                 % propia sola / propia + OneWeb (los dos Ku reales)
alphaStar = 0.5;

fprintf('Operadores (todos en la MISMA banda, %g GHz, sin coordinacion):\n', cfg.radio.freq_GHz);
for i = 1:numel(OPS)
    fprintf('  %d) %-13s h=%4d km  inc=%5.1f deg  %-5s  T=%4d  P=%3d%s\n', i, OPS(i).name, ...
        OPS(i).h, OPS(i).inc, OPS(i).pattern, OPS(i).T, OPS(i).P, ...
        repmat('   <- PROPIO (sirve al usuario)', 1, i==1));
end
fprintf('\nNivel fino: dt = %g s (los eventos in-line duran segundos)\n', dt_fine);
fprintf('Nivel KPI : dt = %g s (INFRAMUESTREA los eventos -> cota inferior)\n\n', cfg.time.dt);

%% ------------------------------------------------------------------------
%  2. TAREA 2: ESTADISTICA DE EVENTOS IN-LINE (nivel fino, M = 1)
%  ------------------------------------------------------------------------
fprintf('=======================================================================\n');
fprintf('  TAREA 2: EVENTOS IN-LINE (dt = %g s, punto de referencia)\n', dt_fine);
fprintf('=======================================================================\n\n');

fprintf('Apuntamiento del haz AJENO: cfg.interf.satPointing = ''%s'' (nominal)\n\n', ...
    cfg.interf.satPointing);

INL = cell(1, numel(nOpsList));
for q = 1:numel(nOpsList)
    nOps = nOpsList(q);
    fprintf('--- %d operador(es): %s ---\n', nOps, strjoin({OPS(1:nOps).name}, ' + '));
    INL{q} = inline_events(cfg, OPS(1:nOps), dt_fine, psiTh);
    fprintf('\n');
end

%% --- COTA SUPERIOR: haz ajeno apuntado AL USUARIO ('boresight') ------------
%  El caso nominal supone que el satelite ajeno apunta su haz a SU PROPIO NADIR, asi
%  que aunque este in-line con nuestra linea de vista nos ilumina con un lobulo
%  lateral. La COTA SUPERIOR realista es que ademas este sirviendo a un usuario
%  cercano al nuestro, es decir que nos apunte con la ganancia de PICO. Es el caso
%  que de verdad preocupa en coordinacion NGSO-NGSO (limites EPFD).
fprintf('=== COTA SUPERIOR: haz ajeno en BORESIGHT (apuntando al usuario) ===\n');
cfgB = cfg;  cfgB.interf.satPointing = 'boresight';
INL_bore = cell(1, numel(nOpsList));
for q = 1:numel(nOpsList)
    if nOpsList(q) == 1, INL_bore{q} = INL{q}; continue; end
    fprintf('--- %d operador(es), BORESIGHT ---\n', nOpsList(q));
    INL_bore{q} = inline_events(cfgB, OPS(1:nOpsList(q)), dt_fine, psiTh);
    fprintf('\n');
end

fprintf('--- NOMINAL (nadir) vs COTA (boresight): caida de SINR en evento psi<%g deg ---\n', psiTh(2));
fprintf('%6s %16s %16s %14s %14s\n','nOps','dSINR nadir[dB]','dSINR bore[dB]','I/N nadir','I/N bore');
for q = 1:numel(nOpsList)
    if nOpsList(q) == 1, continue; end
    fprintf('%6d %16.3f %16.3f %14.2f %14.2f\n', nOpsList(q), ...
        INL{q}.dSINR_event, INL_bore{q}.dSINR_event, ...
        INL{q}.IoN_event,   INL_bore{q}.IoN_event);
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  3. TAREAS 3-4: IMPACTO EN KPIs Y VIABILIDAD (nivel pipeline)
%  ------------------------------------------------------------------------
fprintf('\n=======================================================================\n');
fprintf('  TAREAS 3-4: IMPACTO EN KPIs (pipeline, dt = %g s)\n', cfg.time.dt);
fprintf('=======================================================================\n');
fprintf('[AVISO] dt = %g s INFRAMUESTREA los eventos in-line (duran ~s): lo que sigue\n', cfg.time.dt);
fprintf('        es COTA INFERIOR del impacto. La cola real la da la TAREA 2.\n\n');

cases = { ...
    struct('scheme','reuse1','Delta',1,'alpha',NaN,      'adaptive',false,'name','reuse1'), ...
    struct('scheme','reuseD','Delta',3,'alpha',NaN,      'adaptive',false,'name','reuseD(3)'), ...
    struct('scheme','reuseD','Delta',4,'alpha',NaN,      'adaptive',false,'name','reuseD(4)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',alphaStar,'adaptive',false,'name','ffr(3,a*)'), ...
    struct('scheme','ffr',   'Delta',4,'alpha',alphaStar,'adaptive',false,'name','ffr(4,a*)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',alphaStar,'adaptive',true, 'name','ffr-adapt(3)') };
nC = numel(cases);
nm = cellfun(@(c) c.name, cases, 'UniformOutput', false);

pts = cell(1, numel(nOpsList));
for q = 1:numel(nOpsList)
    pts{q} = struct('constellations', OPS(1:nOpsList(q)), 'cases', {cases}, ...
                    'label', sprintf('%dop', nOpsList(q)));
end

nCores = max(feature('numcores'), 1);
fprintf('--------- PLAN DE MEMORIA ---------\n');
MEM = estimate_sweep_memory(cfg, pts, nCores, struct('verbose',true));
if ~isnan(MEM.timeBlock_reco)
    fprintf('\n[AVISO] Se reduce cfg.compute.timeBlock de %d a %d.\n', ...
        cfg.compute.timeBlock, MEM.timeBlock_reco);
    cfg.compute.timeBlock = MEM.timeBlock_reco;
    MEM = estimate_sweep_memory(cfg, pts, nCores, struct('verbose',false));
end
if isnan(MEM.nFit) || MEM.nFit <= 1, nWorkers = 0; else, nWorkers = min(MEM.nFit, nCores); end
fprintf('\n[plan] pico/worker %.2f GB | caben %d -> %d worker(s)\n\n', ...
    MEM.peak_worker_GB, MEM.nFit, max(nWorkers,1));

[res, INFO] = run_sweep_points(cfg, pts, nWorkers, struct('verbose',true));

%% --- Tablas de KPIs -------------------------------------------------------
fprintf('\n\n--------- KPIs vs NUMERO DE OPERADORES ---------\n');
for c = 1:nC
    fprintf('\n--- %s ---\n', nm{c});
    fprintf('%6s %8s %11s %11s %10s %10s %10s %9s %s\n', ...
        'nOps','N sat','SINRe_p1','SINRe_p5','cob.calid','Rp5_bor','R_agg','penaliz','veredicto');
    for q = 1:numel(nOpsList)
        K = res{q}.K{c};  v = K.viab;
        fprintf('%6d %8d %11.3f %11.3f %10.2f %10.2f %10.2f %9.2f %s\n', ...
            nOpsList(q), res{q}.N, v.SINR_edge_p1, v.SINR_edge_p5, v.coverage, ...
            K.edge.R_p5_Mbps, K.all.R_agg_Mbps/1e3, K.penalty_mean, v.verdict);
    end
end

fprintf('\n--------- CAMBIA ALGUN VEREDICTO al anadir operadores? ---------\n');
anyChange = false;
for c = 1:nC
    v1 = res{1}.K{c}.viab.verdict;
    for q = 2:numel(nOpsList)
        vq = res{q}.K{c}.viab.verdict;
        if ~strcmp(v1, vq)
            fprintf('  %s: %s (1 op) -> %s (%d ops)  ** CAMBIA **\n', nm{c}, v1, vq, nOpsList(q));
            anyChange = true;
        end
    end
end
if ~anyChange
    fprintf('  NINGUNO. Los %d esquemas conservan su veredicto con %s operadores.\n', ...
        nC, strjoin(arrayfun(@(n) sprintf('%d',n), nOpsList, 'UniformOutput', false), ' y '));
    fprintf('  -> la inter-constelacion ENGORDA LA COLA pero no mueve la decision de viabilidad\n');
    fprintf('     (al nivel de muestreo del pipeline; ver la cota de la TAREA 2).\n');
end

fprintf('\n--------- Degradacion por operador anadido (esquema reuse1 y ffr(3)) ---------\n');
fprintf('%6s %12s %12s %12s %12s\n','nOps','r1 p1','r1 p5','ffr3 p1','ffr3 p5');
for q = 1:numel(nOpsList)
    fprintf('%6d %12.3f %12.3f %12.3f %12.3f\n', nOpsList(q), ...
        res{q}.K{1}.viab.SINR_edge_p1, res{q}.K{1}.viab.SINR_edge_p5, ...
        res{q}.K{4}.viab.SINR_edge_p1, res{q}.K{4}.viab.SINR_edge_p5);
end

%% ------------------------------------------------------------------------
%  4. TAREA 4: ESCALADO CON EL NUMERO DE OPERADORES
%  ------------------------------------------------------------------------
fprintf('\n\n=======================================================================\n');
fprintf('  TAREA 4: ESCALADO CON EL NUMERO DE OPERADORES\n');
fprintf('=======================================================================\n\n');
fprintf('%6s %10s %12s', 'nOps','N ajenos','psi_min med');
for x = psiTh, fprintf('   P(psi<%gdeg)', x); end
fprintf(' %12s %12s\n', 'dur med[s]', 'SINR p1');
for q = 1:numel(nOpsList)
    E = INL{q};
    fprintf('%6d %10d %12.3f', nOpsList(q), E.nForeign, E.psi_min_mean);
    for j = 1:numel(psiTh), fprintf(' %13.5f', E.P_event(j)); end
    fprintf(' %12.2f %12.3f\n', E.dur_mean(2), res{q}.K{1}.viab.SINR_edge_p1);
end
fprintf('\n(dur med = duracion media de los eventos con psi < %g deg)\n', psiTh(2));

% NO SE AJUSTA ELASTICIDAD. Con DOS puntos (1 y 2 operadores) una regresion
% log-log tiene cero grados de libertad: pasaria exactamente por los dos puntos y su
% pendiente no seria un resultado, seria una tautologia. La version anterior de este
% experimento la reportaba con TRES puntos, pero el tercero era Kuiper, que opera en
% Ka y no es co-canal en Ku (ver cabecera): esa cifra queda RETIRADA. Lo que SI se
% mide, y con dato real, es el contraste de ARQUITECTURA que viene ahora.

%% --- ORTOGONALIDAD GEOMETRICA: por que OneWeb estorba poco ------------------
%  El hallazgo NO necesita un tercer operador. Se sostiene sobre el contraste entre
%  las DOS arquitecturas reales: OneWeb interfiere poco no porque tenga pocos
%  satelites (720 no es poco) sino porque su geometria es CASI ORTOGONAL a la nuestra.
%  Se cuantifica con el mecanismo dominante, que es su mascara de servicio: OneWeb
%  declara 55 deg y nosotros 25, luego la mayor parte de su constelacion visible sobre
%  nuestro punto NO esta radiando hacia el.
fprintf('\n=======================================================================\n');
fprintf('  ORTOGONALIDAD GEOMETRICA (por que OneWeb estorba poco)\n');
fprintf('=======================================================================\n');
foreign = OPS(2:end);
for f = 1:numel(foreign)
    op = foreign(f);
    visOwn = mean_visible(cfg, op, op.minElev);          % con SU mascara real
    visGlb = mean_visible(cfg, op, cfg.geom.minElev);    % con la mascara global (25)
    fprintf(['  %-12s  T=%d  h=%d km  inc=%.1f  mascara propia %g deg\n' ...
             '                visibles/instante: %.2f con SU mascara  vs  %.2f con %g deg\n' ...
             '                -> la mascara real elimina el %.1f%% de su visibilidad\n'], ...
        op.name, op.T, op.h, op.inc, op.minElev, visOwn, visGlb, cfg.geom.minElev, ...
        100*(1 - visOwn/max(visGlb,eps)));
end
fprintf(['\n  LECTURA: la diversidad de ARQUITECTURA (altitud 1200 vs 550 km, inclinacion\n' ...
         '  87.9 polar vs 53 delta, mascara 55 vs 25 deg) es lo que protege, no la\n' ...
         '  distancia ni el numero de satelites. Un operador Ku con la MISMA altitud y\n' ...
         '  la MISMA mascara que la propia agravaria la cola por este mismo mecanismo,\n' ...
         '  leido al reves. Eso es una PREDICCION CUALITATIVA del modelo, no una cifra\n' ...
         '  simulada: no hay hoy un segundo operador Ku con esa arquitectura.\n']);

%% ------------------------------------------------------------------------
%  5. FIGURAS
%  ------------------------------------------------------------------------
if ~exist('figs_e4','dir'), mkdir('figs_e4'); end
make_e4_figures(INL, nOpsList, psiTh, res, cases, cfg);

%% ------------------------------------------------------------------------
%  6. Guardado
%  ------------------------------------------------------------------------
save('interconstellation_results.mat', 'INL', 'INL_bore', 'res', 'INFO', 'OPS', ...
     'nOpsList', 'psiTh', 'cases', 'nm', 'cfg', 'dt_fine', 'MEM', 'alphaStar', '-v7.3');
fprintf('\nResultados en interconstellation_results.mat | figuras en figs_e4/\n');

%% ========================================================================
%  FUNCIONES LOCALES
%  ========================================================================
function v = mean_visible(cfgBase, op, minElevDeg)
%MEAN_VISIBLE  Numero medio de satelites de UNA constelacion visibles por instante.
%   Solo geometria (no radioenlace): construye la constelacion `op`, la propaga sobre
%   la ventana del experimento y cuenta cuantos superan `minElevDeg` desde el punto de
%   referencia, promediando en el tiempo. Se usa para cuantificar cuanta de la
%   constelacion ajena esta REALMENTE radiando hacia nosotros con su mascara real
%   frente a la mascara global. Barato: M = 1.
cfg = cfgBase;
cfg.constellations = op;
cfg.constellations(1).minElev = minElevDeg;
cfg.ground.mode = 'point';
cfg.time.dt     = 60;                 % basta para una MEDIA temporal

tvec  = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
sats  = build_constellation(cfg);
users = build_user_grid(cfg);

nVis = nan(1, numel(tvec));
blk  = max(1, floor(2e6 / max(sats.N,1)));
for b = 1:ceil(numel(tvec)/blk)
    kk = ((b-1)*blk+1) : min(b*blk, numel(tvec));
    tv = tvec(kk);
    R  = eci2ecef(cfg, propagate(cfg, sats, tv), tv);
    G  = compute_geometry(cfg, users, R, sats.minElev);
    nVis(kk) = squeeze(sum(G.vis(1,:,:), 2)).';
    clear R G;
end
v = mean(nVis, 'omitnan');
end

% -------------------------------------------------------------------------
function E = inline_events(cfgBase, ops, dt, psiTh)
%INLINE_EVENTS  Estadistica de eventos in-line en el punto de referencia.
%   Nivel FINO: M = 1 y dt de segundos. Replica EXACTAMENTE las formulas de
%   compute_interference (separacion angular psi por la esferica, patron del haz
%   ajeno segun cfg.interf.satPointing con suelo de lobulos, terminal ITU-R S.1428
%   via term_gain_dB), pero SEPARANDO la contribucion de cada constelacion, que es lo
%   que compute_interference no distingue.
%
%   Devuelve, sobre los instantes CON COBERTURA del operador propio:
%     .psi_min        [1 x Nt] separacion angular al satelite AJENO mas alineado
%     .P_event        P(psi_min < X) para cada umbral X
%     .dur_mean/.dur_max  duracion media/maxima de los eventos [s]
%     .n_events       numero de eventos por umbral
%     .SINR_with/.SINR_without   SINR con y sin la interferencia ajena [dB]
%     .dSINR_event / .dSINR_quiet  caida media de SINR dentro/fuera de evento [dB]
%     .IoN_foreign    I_ajena / N [dB] por instante

cfg = cfgBase;
cfg.constellations = ops;
cfg.ground.mode = 'point';                 % M = 1: el coste es N*Nt
cfg.time.dt = dt;

tvec = cfg.time.t0 : dt : cfg.time.t0 + cfg.time.duration;
Nt   = numel(tvec);
sats = build_constellation(cfg);
users = build_user_grid(cfg);
N    = sats.N;
own  = sats.cid == 1;
E.nForeign = sum(~own);

Re         = cfg.const.Re;
EIRP_const = cfg.radio.EIRPdensity_dBWMHz + 10*log10(cfg.radio.B_MHz);
Grx_max    = cfg.radio.Grx_max_dBi;
N_dBW      = -228.6 + 10*log10(cfg.radio.Tsys_K) + 10*log10(radio_band_Hz(cfg));
N_lin      = 10^(N_dBW/10);

psi_min = nan(1,Nt);  Iown = zeros(1,Nt);  Ifor = zeros(1,Nt);
Clin    = nan(1,Nt);  covOK = false(1,Nt);

% Troceado temporal: con dt = 1 s y N ~ miles, los tensores N*Nt no caben de una pieza.
blk = max(1, floor(2e6 / max(N,1)));       % ~2e6 elementos por bloque
nBlocks = ceil(Nt/blk);
for b = 1:nBlocks
    kk = ((b-1)*blk+1) : min(b*blk, Nt);
    tv = tvec(kk);
    R_ecef = eci2ecef(cfg, propagate(cfg, sats, tv), tv);
    % Mascara POR CONSTELACION: OneWeb no opera por debajo de 55 deg, luego sus
    % satelites por debajo de esa elevacion NO deben contar como interferentes.
    G = compute_geometry(cfg, users, R_ecef, sats.minElev);
    % Servidor SOLO entre los propios. Se usa el MOTOR COMPARTIDO con su satMask, no
    % una copia local: asi el nivel FINO y el nivel KPI aplican EXACTAMENTE la misma
    % regla de asociacion, incluida cfg.geom.serving ('maxElev' | 'minRange').
    % (auditoria, hallazgo B12: la copia local `associate_serving_own` hardcodeaba
    % maxElev e ignoraba cfg.geom.serving; si alguien cambiase el criterio, los dos
    % niveles habrian divergido en silencio -- que es exactamente la situacion que
    % destapo el bug de servidor ajeno de E4. Verificado antes de borrarla:
    % S.idx identico y max|S.el dif| = 0 sobre 3388 satelites y 61 instantes.)
    S = associate_serving(cfg, G, own);
    FSPL = freespace_loss(G.range, cfg.radio.freq_GHz);
    Latm = atm_loss_dB(cfg, G.el);

    for j = 1:numel(kk)
        k = kk(j);
        A = S.idx(j);
        if isnan(A), continue; end
        covOK(k) = true;

        C_dBW  = EIRP_const + Grx_max - FSPL(1,A,j) - Latm(1,A,j);
        Clin(k) = 10^(C_dBW/10);

        vis = squeeze(G.vis(1,:,j)).';
        vis(A) = false;
        if ~any(vis), continue; end

        el_A = G.el(1,A,j);  az_A = G.az(1,A,j);
        idx  = find(vis);
        el_B = reshape(G.el(1,idx,j), 1, []);
        az_B = reshape(G.az(1,idx,j), 1, []);
        h_B  = reshape(sats.h(idx),   1, []);
        F_B  = reshape(FSPL(1,idx,j), 1, []);
        L_B  = reshape(Latm(1,idx,j), 1, []);

        % Separacion angular usuario->servidor vs usuario->interferente (esferica)
        carg = sind(el_A).*sind(el_B) + cosd(el_A).*cosd(el_B).*cosd(az_A - az_B);
        psi  = acosd(max(-1, min(1, carg)));

        % Patron del haz del interferente (mismo criterio que compute_interference)
        switch lower(cfg.interf.satPointing)
            case 'boresight'
                Gsat_rel = zeros(1, numel(idx));
            otherwise
                eta = asind( (Re ./ (Re + h_B)) .* cosd(el_B) );
                Gsat_rel = beam_gain_dB(eta, cfg.radio.beamwidth3dB_deg, ...
                                        cfg.radio.sidelobe_floor_dB);
        end
        Grx_B = term_gain_dB(cfg, psi);                    % ITU-R S.1428
        I_B   = 10.^((EIRP_const + Gsat_rel + Grx_B - F_B - L_B)/10);

        isFor = ~own(idx).';
        Ifor(k) = sum(I_B(isFor));
        Iown(k) = sum(I_B(~isFor));
        if any(isFor), psi_min(k) = min(psi(isFor)); end
    end
    clear R_ecef G S FSPL Latm;
end

% --- Estadistica ---------------------------------------------------------
E.tvec     = tvec;
E.psi_min  = psi_min;
E.cov      = covOK;
E.IoN_foreign = 10*log10(max(Ifor,realmin)/N_lin);
E.SINR_with    = 10*log10(Clin ./ (N_lin + Iown + Ifor));
E.SINR_without = 10*log10(Clin ./ (N_lin + Iown));
E.dSINR        = E.SINR_without - E.SINR_with;         % caida por la ajena [dB]

valid = covOK & ~isnan(psi_min);
E.n_valid      = sum(valid);
E.P_event  = zeros(1,numel(psiTh));
E.dur_mean = zeros(1,numel(psiTh));  E.dur_max = zeros(1,numel(psiTh));
E.n_events = zeros(1,numel(psiTh));

% Caso de referencia SIN operadores ajenos: no hay psi que medir. Se devuelve la
% estructura completa con NaN/0 para que las tablas comparativas no se rompan.
if E.nForeign == 0 || E.n_valid == 0
    E.psi_min_mean = NaN;  E.psi_min_min = NaN;
    E.dSINR_quiet = 0;  E.dSINR_event = 0;
    E.IoN_event = -Inf;  E.IoN_quiet = -Inf;
    if any(covOK)
        E.SINR_p1_fine = prctile(E.SINR_with(covOK), 1);
    else
        E.SINR_p1_fine = NaN;
    end
    E.SINR_p1_fine_noForeign = E.SINR_p1_fine;
    fprintf('  N ajenos = 0 -> CASO DE REFERENCIA (solo la constelacion propia).\n');
    fprintf('  SINR p1 (nivel fino) = %.3f dB. No hay eventos in-line que medir.\n', E.SINR_p1_fine);
    return;
end

E.psi_min_mean = mean(psi_min(valid));
E.psi_min_min  = min(psi_min(valid));
fprintf('  N ajenos = %d | instantes validos = %d | psi_min medio = %.2f deg (min %.3f)\n', ...
    E.nForeign, E.n_valid, E.psi_min_mean, E.psi_min_min);
fprintf('  %8s %12s %10s %12s %12s %14s\n','psi<X','P(evento)','n eventos','dur med[s]','dur max[s]','dSINR med[dB]');
for j = 1:numel(psiTh)
    ev = valid & (psi_min < psiTh(j));
    E.P_event(j) = sum(ev) / max(sum(valid),1);
    [nEv, dMean, dMax] = run_lengths(ev, dt);
    E.n_events(j) = nEv;  E.dur_mean(j) = dMean;  E.dur_max(j) = dMax;
    if any(ev), dS = mean(E.dSINR(ev)); else, dS = 0; end
    fprintf('  %8.1f %12.5f %10d %12.2f %12.2f %14.3f\n', ...
        psiTh(j), E.P_event(j), nEv, dMean, dMax, dS);
end

quiet = valid & (psi_min >= psiTh(end));
E.dSINR_quiet = mean(E.dSINR(quiet));
E.dSINR_event = mean(E.dSINR(valid & psi_min < psiTh(2)));
E.IoN_event   = mean(E.IoN_foreign(valid & psi_min < psiTh(2)));
E.IoN_quiet   = mean(E.IoN_foreign(quiet));
E.SINR_p1_fine  = prctile(E.SINR_with(valid), 1);
E.SINR_p1_fine_noForeign = prctile(E.SINR_without(valid), 1);
fprintf('  Caida media de SINR:  DURANTE evento (psi<%g) = %.3f dB | fuera (psi>=%g) = %.4f dB\n', ...
    psiTh(2), E.dSINR_event, psiTh(end), E.dSINR_quiet);
fprintf('  I_ajena/N:            DURANTE evento = %.2f dB | fuera = %.2f dB\n', ...
    E.IoN_event, E.IoN_quiet);
fprintf('  SINR p1 (nivel fino): con ajenos = %.3f dB | sin ajenos = %.3f dB  (cola: %.3f dB)\n', ...
    E.SINR_p1_fine, E.SINR_p1_fine_noForeign, E.SINR_p1_fine_noForeign - E.SINR_p1_fine);
end

% -------------------------------------------------------------------------
% NOTA (auditoria, hallazgo B12): aqui vivia `associate_serving_own`, una copia local
% que elegia el servidor entre los satelites propios. Se ha BORRADO: `associate_serving`
% acepta un `satMask` opcional que hace exactamente lo mismo, y ademas respeta
% cfg.geom.serving en vez de hardcodear 'maxElev'. Mantener la copia era el patron
% "dos reglas para lo mismo" que en este proyecto ya causo los bugs de atm_key_of y
% del servidor de otra constelacion (este mismo experimento). Equivalencia verificada
% antes de borrarla: S.idx identico y max|S.el dif| = 0 sobre 3388 satelites x 61
% instantes con los 3 operadores.

% -------------------------------------------------------------------------
function [nEv, dMean, dMax] = run_lengths(mask, dt)
%RUN_LENGTHS  Numero y duracion de las rachas TRUE consecutivas de `mask`.
m = [false, mask(:).', false];
d = diff(m);
ini = find(d == 1);  fin = find(d == -1);
len = (fin - ini) * dt;
nEv = numel(len);
if nEv == 0, dMean = 0; dMax = 0; else, dMean = mean(len); dMax = max(len); end
end

% -------------------------------------------------------------------------
function make_e4_figures(INL, nOpsList, psiTh, res, cases, cfg)
%MAKE_E4_FIGURES  Figuras de E4.
nQ = numel(nOpsList);  nC = numel(cases);
nm = cellfun(@(c) c.name, cases, 'UniformOutput', false);
col = lines(max(nC,nQ));

% (a) psi_min vs tiempo junto a la SINR: los picos deben coincidir
E = INL{end};                                  % el caso con mas operadores
v = E.cov & ~isnan(E.psi_min);
figure('Name','E4(a) Eventos in-line en el tiempo','Position',[60 60 1100 480]);
yyaxis left;
plot(E.tvec(v)/60, E.psi_min(v), 'LineWidth', 0.8); hold on;
yline(psiTh(2), '--', sprintf('%g deg', psiTh(2)), 'LineWidth', 1.2);
ylabel('\psi_{min} al satelite ajeno mas alineado [deg]'); set(gca,'YScale','log');
yyaxis right;
plot(E.tvec(v)/60, E.SINR_with(v), 'LineWidth', 1.0);
ylabel('SINR [dB]');
xlabel('tiempo [min]'); grid on;
title(sprintf('E4(a): eventos in-line (%d operadores) - los minimos de \\psi coinciden con los picos de interferencia', nOpsList(end)));
save_fig(gcf, 'figs_e4', 'e4_a_inline_timeline');

% (b) CDF de psi_min por numero de operadores
figure('Name','E4(b) CDF de psi_min','Position',[60 60 1100 460]);
subplot(1,2,1); hold on; grid on;
for q = 1:nQ
    Eq = INL{q};  vq = Eq.cov & ~isnan(Eq.psi_min);
    if ~any(vq), continue; end
    x = sort(Eq.psi_min(vq));
    plot(x, (1:numel(x))/numel(x), 'LineWidth', 1.6, 'Color', col(q,:));
end
set(gca,'XScale','log');
for j = 1:numel(psiTh), xline(psiTh(j), ':', sprintf('%g', psiTh(j))); end
xlabel('\psi_{min} [deg]'); ylabel('F(\psi_{min})'); title('CDF de la separacion al ajeno mas alineado');
legend(arrayfun(@(n) sprintf('%d op', n), nOpsList, 'UniformOutput', false), 'Location','southeast');

subplot(1,2,2); hold on; grid on;
for j = 1:numel(psiTh)
    plot(nOpsList, arrayfun(@(q) INL{q}.P_event(j), 1:nQ), '-o', 'LineWidth', 1.6);
end
set(gca,'YScale','log');
xlabel('numero de operadores'); ylabel('P(\psi_{min} < X)');
title('Probabilidad de evento in-line');
legend(arrayfun(@(x) sprintf('X = %g deg', x), psiTh, 'UniformOutput', false), 'Location','southeast');
sgtitle('E4(b): estadistica de alineamientos');
save_fig(gcf, 'figs_e4', 'e4_b_psimin_cdf');

% (c) Impacto en KPIs: p1 y p5 vs numero de operadores
figure('Name','E4(c) KPIs vs operadores','Position',[60 60 1100 460]);
subplot(1,2,1); hold on; grid on;
for c = 1:nC
    plot(nOpsList, arrayfun(@(q) res{q}.K{c}.viab.SINR_edge_p1, 1:nQ), '-o', ...
        'Color', col(c,:), 'LineWidth', 1.6);
end
yline(cfg.viab.gamma_th_dB,'k--'); yline(cfg.viab.gamma_floor_dB,'r--');
xlabel('numero de operadores'); ylabel('SINR_{edge,p1} [dB]');
title('Percentil 1 (la COLA)'); legend(nm,'Location','best','FontSize',7);
subplot(1,2,2); hold on; grid on;
for c = 1:nC
    plot(nOpsList, arrayfun(@(q) res{q}.K{c}.viab.SINR_edge_p5, 1:nQ), '-o', ...
        'Color', col(c,:), 'LineWidth', 1.6);
end
yline(cfg.viab.gamma_th_dB,'k--'); yline(cfg.viab.gamma_floor_dB,'r--');
xlabel('numero de operadores'); ylabel('SINR_{edge,p5} [dB]');
title('Percentil 5 (metrica reina)');
sgtitle('E4(c): impacto en KPIs (COTA INFERIOR: dt del pipeline inframuestrea los eventos)');
save_fig(gcf, 'figs_e4', 'e4_c_kpi_vs_operators');

% (d) Caida de SINR y I/N dentro vs fuera de evento
figure('Name','E4(d) Efecto de los eventos','Position',[60 60 1100 460]);
subplot(1,2,1); hold on; grid on;
Ev = INL{end};  vv = Ev.cov & ~isnan(Ev.psi_min);
scatter(Ev.psi_min(vv), Ev.dSINR(vv), 6, 'filled', 'MarkerFaceAlpha', 0.25);
set(gca,'XScale','log');
xlabel('\psi_{min} [deg]'); ylabel('caida de SINR por la ajena [dB]');
title('La caida se dispara al alinearse');
subplot(1,2,2); hold on; grid on;
scatter(Ev.psi_min(vv), Ev.IoN_foreign(vv), 6, 'filled', 'MarkerFaceAlpha', 0.25);
set(gca,'XScale','log'); yline(0,'k--','I = N');
xlabel('\psi_{min} [deg]'); ylabel('I_{ajena} / N [dB]');
title('Interferencia ajena frente al ruido');
sgtitle('E4(d): mecanismo - el terminal S.1428 deja de discriminar cuando \psi \rightarrow 0');
save_fig(gcf, 'figs_e4', 'e4_d_event_effect');

fprintf('\n4 figuras exportadas a figs_e4/\n');
end
