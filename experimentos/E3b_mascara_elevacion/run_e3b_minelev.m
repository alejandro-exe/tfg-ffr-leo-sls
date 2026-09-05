%% run_e3b_minelev  E3b: COMPRESION angular -- el eje que SI rompe la auto-similitud.
%
%  POR QUE ESTE EJE Y NO EL DE E3
%  -------------------------------
%  E3 (densificacion de haces) NO satura el borde, y la causa es una IDENTIDAD, no
%  una casualidad numerica. Con la separacion normativa 3GPP TR 38.821 Sec. 6.1.1
%
%       s = sqrt(3) * sin(HPBW/2)            [rad]
%
%  y el patron Bessel de 38.811 Sec. 6.4.1 parametrizado por el ancho de haz
%
%       u = 1.61634 * sin(theta_off) / sin(HPBW/2)
%
%  el argumento hacia el vecino CO-CANAL DEL PRIMER ANILLO (theta_off = s) vale
%
%       u = 1.61634 * sin(s)/sin(HPBW/2)  ~  1.61634 * sqrt(3)  =  2.79957
%       G = 20*log10|2*J1(u)/u|           =  -10.67 dB
%
%  INDEPENDIENTE DEL ANCHO DE HAZ. La reticula co-canal es INVARIANTE DE ESCALA: al
%  estrechar el haz, celdas y separacion encogen JUNTAS y el vecino sigue estando a
%  los mismos 0.866 anchos de haz. Por eso densificar multiplica la capacidad sin
%  degradar la C/I (E3: -0.85 dB de SINR_edge_p5 por x36 de densidad). El script
%  VERIFICA esta identidad numericamente de 3 a 0.1 deg antes de empezar (seccion 1).
%
%  MATIZ IMPORTANTE AL CITAR ESA CIFRA. Los -10.67 dB son los del vecino
%  ADYACENTE (theta_off = s), que es el primer co-canal SOLO EN REUSO PLENO, donde
%  todos los haces comparten banda. Con reuso Delta el primer co-canal esta mas lejos
%  y cae mas abajo (limite de angulo pequeno):
%       Delta=3 -> sqrt(3)*s = 1.500 anchos de haz -> -17.954 dB
%       Delta=4 -> 2*s       = 1.732 anchos de haz -> -18.456 dB
%  La PROPIEDAD de invariancia se cumple igual de bien en las TRES distancias
%  (recorridos de 0.010 / 0.014 / 0.028 dB sobre el eje de E3), luego lo que hay que
%  acotar es EL NUMERO, no la propiedad: la conclusion "densificar no satura" vale
%  para los tres esquemas.
%  TRAMPA: el sqrt(3) de u = 1.61634*sqrt(3) procede de la separacion normativa
%  s = sqrt(3)*sin(HPBW/2), NO de la distancia de reuso sqrt(3)*s de Delta=3. La
%  coincidencia es casual y en una defensa oral induce a error: u(s) = 2.79958 frente
%  a u(sqrt(3)*s) = 4.84902.
%
%  La ELEVACION rompe esa invariancia porque NO es un cambio de escala. Las celdas
%  son EARTH-FIXED: su teselado esta fijo en el suelo, y lo que varia es desde donde
%  se mira. A elevacion baja el satelite ve el cluster ESCORZADO: la separacion
%  angular efectiva entre centros de celda vecinos se COMPRIME en la direccion
%  radial (hacia el punto subsatelite) por un factor ~sin(elev), mientras el ancho
%  de haz HPBW sigue siendo el mismo. Resultado: los haces vecinos se solapan mas,
%  u baja, la ganancia co-canal SUBE y la C/I se hunde. Es una compresion GEOMETRICA
%  que ningun reescalado de las celdas compensa.
%
%  Magnitud del efecto ya medida en E5 (desglose por tramo de elevacion):
%  SINR_p5 de borde -11.58 dB en 25-35 deg frente a -4.18 dB en 60-90 deg, es decir
%  ~7 dB, frente a ~1 dB de TODO el eje de densificacion de E3.
%
%  EJE: cfg.geom.minElev = [45 40 35 30 25 20 15 10] deg. Relajar la mascara amplia
%  la cobertura geometrica y la disponibilidad, a costa de admitir geometrias peores.
%  Todo lo demas FIJO (config convergida de FASE A; HPBW=1.5, nRings=5).
%
%  NOTA DE COSTE: cfg.geom.minElev forma parte de la CLAVE de la cache P.618 de
%  atm_loss_dB (freq|lat|lon|minElev|p618_avail), asi que cada punto necesita su
%  PROPIA tabla (~85 s) en vez de compartir una como en E3. run_sweep_points lo
%  detecta y lo avisa; el barrido tarda bastante mas de lo que su tamano sugiere.
%
%  Uso: ejecutar en la raiz del proyecto. Resultados en minelev_results.mat,
%  figuras en figs_e3b/. No modifica la fisica del motor.

clear; clc; close all;

fprintf('=======================================================================\n');
fprintf('  E3b: COMPRESION ANGULAR (barrido de minElev)\n');
fprintf('=======================================================================\n\n');

%% ------------------------------------------------------------------------
%  1. VERIFICACION de la identidad de invariancia de escala (motivacion del eje)
%  ------------------------------------------------------------------------
% MATIZ: la linea que imprime esta cabecera dice "vecino co-canal del 1er anillo", y
% eso SOLO es cierto en REUSO PLENO: ese vecino es el ADYACENTE (theta_off = s). Con
% Delta=3 el primer co-canal esta a sqrt(3)*s (-17.954 dB) y con Delta=4 a 2*s
% (-18.456 dB). El texto impreso no se reescribe a proposito, para no alterar la
% salida ya publicada del experimento.
fprintf('--------- Invariancia de escala de la reticula co-canal ---------\n');
fprintf('Vecino co-canal del 1er anillo: theta_off = s = sqrt(3)*sin(HPBW/2)\n');
fprintf('%10s %12s %12s %12s %12s\n','HPBW[deg]','s[deg]','u','G_1er[dB]','u/2.79957');
uref = 1.61634*sqrt(3);
hpbwChk = [3 2.2 1.5 1.0 0.7 0.5 0.3 0.1];
Gchk = zeros(size(hpbwChk));
for i = 1:numel(hpbwChk)
    hb  = hpbwChk(i);
    sd  = rad2deg(sqrt(3)*sin(deg2rad(hb)/2));
    uu  = 1.61634*sind(sd)/sind(hb/2);
    Gchk(i) = beam_gain_dB(sd, hb, -Inf);
    fprintf('%10.2f %12.4f %12.5f %12.3f %12.6f\n', hb, sd, uu, Gchk(i), uu/uref);
end
fprintf('u teorico = 1.61634*sqrt(3) = %.5f  ->  G = %.3f dB\n', uref, ...
    20*log10(abs(2*besselj(1,uref)/uref)));
fprintf('Dispersion de G en 3..0.1 deg: %.4f dB  -> INVARIANTE (por eso E3 no satura)\n\n', ...
    max(Gchk)-min(Gchk));

%% ------------------------------------------------------------------------
%  2. CONFIGURACION BASE (perfil convergido de FASE A, todo FIJO salvo minElev)
%  ------------------------------------------------------------------------
cfg = config_default();

cfg.constellations(1).T = 1584;
cfg.constellations(1).P = 72;

cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 4;
cfg.time.dt          = 120;
cfg.time.duration    = 7200;

cfg.ffr.tau_mode = 'quantile';
cfg.ffr.tau_q    = 50;

cfg.beams.nRings = 5;
cfg.beams.nBeams = [];

cfg.compute.timeBlock = 20;

minElevList = [45 40 35 30 25 20 15 10];
alphaStar   = 0.5;
nWorkers    = [];

% DOS DENSIDADES (ver seccion 3b, "binding" de la mascara). La mascara de elevacion
% solo controla la geometria del SERVIDOR cuando llega a MORDER, y con
% cfg.geom.serving = 'maxElev' el servidor es SIEMPRE el satelite mas alto: con una
% constelacion densa el mejor satelite nunca baja de ~52 deg y la mascara es
% IRRELEVANTE para el servidor (solo cambia cuantos INTERFERENTES se ven).
%   'denso'    T=1584 -> el barrido especificado; sale DEGENERADO (se documenta).
%   'disperso' T=66   -> perfil de run_ffr_demo/E5, donde la mascara SI muerde y el
%                        eje de elevacion existe de verdad.
% La clave de la cache P.618 (freq|lat|lon|minElev|p618_avail) NO incluye T, asi que
% las dos densidades COMPARTEN tabla y el coste no se duplica.
densList = struct( ...
    'name', {'denso','disperso'}, ...
    'T',    {1584, 66}, ...
    'P',    {72, 6});
nD = numel(densList);

cases = { ...
    struct('scheme','reuse1','Delta',1,'alpha',NaN,      'adaptive',false,'name','reuse1'), ...
    struct('scheme','reuseD','Delta',3,'alpha',NaN,      'adaptive',false,'name','reuseD(3)'), ...
    struct('scheme','reuseD','Delta',4,'alpha',NaN,      'adaptive',false,'name','reuseD(4)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',alphaStar,'adaptive',false,'name','ffr(3,a*)'), ...
    struct('scheme','ffr',   'Delta',4,'alpha',alphaStar,'adaptive',false,'name','ffr(4,a*)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',alphaStar,'adaptive',true, 'name','ffr-adapt(3)') };
nC = numel(cases);
nm = cellfun(@(c) c.name, cases, 'UniformOutput', false);
nE = numel(minElevList);

fprintf('Perfil: T=%d P=%d | rejilla %.0f km paso %.0f km | dt=%.0f s | HPBW=%.2f | nRings=%d\n', ...
    cfg.constellations(1).T, cfg.constellations(1).P, cfg.ground.radius_km, ...
    cfg.ground.step_km, cfg.time.dt, cfg.radio.beamwidth3dB_deg, cfg.beams.nRings);
fprintf('Eje: minElev = [%s] deg | alpha* = %.2f | %d esquemas\n\n', ...
    num2str(minElevList), alphaStar, nC);

%% ------------------------------------------------------------------------
%  3. GUARDA (invariante en este eje) y COMPRESION ANGULAR ANALITICA
%  ------------------------------------------------------------------------
BL0 = build_beam_layout(cfg);
dCell   = sqrt(sum(BL0.cell_km.^2, 2));
dOutMin = min(dCell(BL0.ring == cfg.beams.nRings));
guard   = (dOutMin - cfg.ground.radius_km) / BL0.spacing_km;

fprintf('--------- GUARDA ---------\n');
fprintf('El layout de haces NO depende de minElev (build_beam_layout solo usa HPBW,\n');
fprintf('nRings y la altitud), luego la guarda es la MISMA en los %d puntos:\n', nE);
fprintf('  s = %.4f deg = %.2f km | anillo externo a %.1f km | rejilla %.0f km -> guarda %.2f anillos\n', ...
    BL0.spacing_deg, BL0.spacing_km, dOutMin, cfg.ground.radius_km, guard);
if guard >= 1
    fprintf('  -> SUFICIENTE (>= 1 anillo) en TODOS los puntos: ningun punto sesgado por guarda.\n');
else
    fprintf('  -> INSUFICIENTE: TODOS los puntos quedarian sesgados.\n');
end
fprintf(['  Se comprueba ademas la senal OPERATIVA por punto (max theta del usuario a su\n' ...
         '  centro de haz frente a s): si la rejilla sobresaliera, ffr_context avisaria.\n\n']);

% --- Compresion angular ANALITICA: separacion efectiva entre celdas vecinas ------
% Dos centros de celda adyacentes separados s_km sobre el suelo, vistos desde el
% satelite a elevacion el. Se evalua en la direccion RADIAL (la que se escorza) y en
% la TRANSVERSAL (que no se escorza), para ver la anisotropia.
elGrid = 5:1:90;
[sepRad, sepTan] = deal(zeros(size(elGrid)));
for i = 1:numel(elGrid)
    [sepRad(i), sepTan(i)] = cell_sep_eff(elGrid(i), BL0.spacing_km, ...
                                          cfg.constellations(1).h, cfg.const.Re);
end
sepNom = BL0.spacing_deg;

fprintf('--------- COMPRESION ANGULAR (analitica) ---------\n');
fprintf('Separacion angular efectiva entre centros de celda vecinos vista desde el satelite\n');
fprintf('(nominal a nadir: s = %.4f deg)\n', sepNom);
fprintf('%10s %14s %14s %12s %14s\n','elev[deg]','sep radial','sep transv','rad/nom','u_1er_anillo');
for el = [90 60 45 35 25 20 15 10]
    [sr, st] = cell_sep_eff(el, BL0.spacing_km, cfg.constellations(1).h, cfg.const.Re);
    uu = 1.61634*sind(sr)/sind(cfg.radio.beamwidth3dB_deg/2);
    fprintf('%10d %14.4f %14.4f %12.3f %14.4f\n', el, sr, st, sr/sepNom, uu);
end
fprintf(['-> Al bajar la elevacion la separacion RADIAL se comprime (~sin(elev)) mientras\n' ...
         '   HPBW no cambia: u baja, la ganancia del vecino co-canal SUBE y la C/I cae.\n' ...
         '   La transversal apenas cambia: la compresion es ANISOTROPA.\n\n']);

%% ------------------------------------------------------------------------
%  3b. DIAGNOSTICO DE "BINDING": la mascara controla la geometria del servidor?
%  ------------------------------------------------------------------------
%  Solo geometria (sin radioenlace ni P.618): barato. Es el diagnostico que decide
%  si el eje de este experimento existe en cada densidad.
fprintf('--------- MUERDE LA MASCARA? (diagnostico de binding) ---------\n');
fprintf(['cfg.geom.serving = ''%s'': el servidor es el satelite MAS ALTO visible.\n' ...
         'Si el mejor satelite nunca baja de la mascara, esta NO afecta al servidor.\n\n'], ...
    cfg.geom.serving);
cfgPt = cfg;  cfgPt.ground.mode = 'point';
tvD   = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
uD    = build_user_grid(cfgPt);
BIND  = struct();
for d = 1:nD
    cD = cfgPt;
    cD.constellations(1).T = densList(d).T;
    cD.constellations(1).P = densList(d).P;
    satD = build_constellation(cD);
    RD   = eci2ecef(cD, propagate(cD, satD, tvD), tvD);
    fprintf('  T=%d:\n', densList(d).T);
    fprintf('  %9s %10s %11s %10s %10s %10s\n','minElev','cobertura','elev_media','elev_p5','elev_min','nVis_med');
    cv = zeros(1,nE); em = zeros(1,nE); e5 = zeros(1,nE); emin = zeros(1,nE); nv = zeros(1,nE);
    for i = 1:nE
        cI = cD;  cI.geom.minElev = minElevList(i);
        GI = compute_geometry(cI, uD, RD);
        SI = associate_serving(cI, GI);
        ok = ~isnan(SI.el);  e = SI.el(ok);
        cv(i)=mean(ok(:)); em(i)=mean(e); e5(i)=prctile(e,5); emin(i)=min(e); nv(i)=mean(SI.nVis(:));
        fprintf('  %9d %10.3f %11.2f %10.2f %10.2f %10.2f\n', ...
            minElevList(i), cv(i), em(i), e5(i), emin(i), nv(i));
    end
    BIND(d).name = densList(d).name;  BIND(d).T = densList(d).T;
    BIND(d).cov = cv; BIND(d).elev_mean = em; BIND(d).elev_p5 = e5;
    BIND(d).elev_min = emin; BIND(d).nVis = nv;
    % La mascara MUERDE si mover el limite cambia la elevacion del servidor
    BIND(d).binds = (max(em) - min(em)) > 0.5;
    if BIND(d).binds
        fprintf('  -> La mascara MUERDE (elev media varia %.2f deg): el eje EXISTE.\n\n', max(em)-min(em));
    else
        fprintf(['  -> La mascara NO MUERDE (elev media varia %.2f deg): el mejor satelite\n' ...
                 '     nunca baja de %.1f deg, luego el barrido sale DEGENERADO en esta\n' ...
                 '     densidad. Solo cambia nVis (%.1f -> %.1f), es decir el numero de\n' ...
                 '     INTERFERENTES inter-satelite, que la FASE A ya probo despreciable.\n\n'], ...
            max(em)-min(em), min(emin), min(nv), max(nv));
    end
end

%% ------------------------------------------------------------------------
%  4. PUNTOS Y PLAN DE MEMORIA
%  ------------------------------------------------------------------------
% nD densidades x nE mascaras. T va en pt.T/pt.P (run_one_density lo aplica DESPUES
% del cfg_over), y minElev en cfg_over.
pts = cell(1, nD*nE);
ptMeta = zeros(nD*nE, 2);         % [densidad, mascara]
k = 0;
for d = 1:nD
    for i = 1:nE
        k = k + 1;
        pts{k} = struct( ...
            'cfg_over', struct('geom', struct('minElev', minElevList(i))), ...
            'T', densList(d).T, 'P', densList(d).P, ...
            'cases',    {cases}, ...
            'label',    sprintf('%s-el%02d', densList(d).name, minElevList(i)));
        ptMeta(k,:) = [d, i];
    end
end

nCores = max(feature('numcores'), 1);
fprintf('--------- PLAN DE MEMORIA ---------\n');
MEM = estimate_sweep_memory(cfg, pts, nCores, struct('verbose',true));
if ~isnan(MEM.timeBlock_reco)
    fprintf('\n[AVISO] Se reduce cfg.compute.timeBlock de %d a %d para que quepan los workers.\n', ...
        cfg.compute.timeBlock, MEM.timeBlock_reco);
    cfg.compute.timeBlock = MEM.timeBlock_reco;
    MEM = estimate_sweep_memory(cfg, pts, nCores, struct('verbose',false));
end
if isempty(nWorkers)
    if isnan(MEM.nFit) || MEM.nFit <= 1, nWorkers = 0; else, nWorkers = min(MEM.nFit, nCores); end
end
fprintf('\n[plan] pico/worker %.2f GB | caben %d -> %d worker(s)\n', ...
    MEM.peak_worker_GB, MEM.nFit, max(nWorkers,1));
fprintf(['[coste] minElev esta en la CLAVE de la cache P.618: los %d puntos necesitan %d\n' ...
         '        tablas DISTINTAS (~85 s cada una), pero T NO esta en la clave, asi que las\n' ...
         '        %d densidades las COMPARTEN. run_sweep_points lo verifica abajo.\n\n'], ...
         numel(pts), nE, nD);

%% ------------------------------------------------------------------------
%  5. BARRIDO
%  ------------------------------------------------------------------------
[res, INFO] = run_sweep_points(cfg, pts, nWorkers, struct('verbose',true));

%% ------------------------------------------------------------------------
%  6. TABLAS
%  ------------------------------------------------------------------------
fprintf('\n\n=======================================================================\n');
fprintf('  TAREA 1/3: KPIs POR ESQUEMA\n');
fprintf('=======================================================================\n');

% Indice: R{d}{i} = resultado de la densidad d, mascara i
R = cell(1,nD);
for d = 1:nD
    R{d} = cell(1,nE);
    for k = 1:numel(pts)
        if ptMeta(k,1) == d, R{d}{ptMeta(k,2)} = res{k}; end
    end
end

for d = 1:nD
    fprintf('\n\n#################  DENSIDAD %s (T = %d)  #################\n', ...
        upper(densList(d).name), densList(d).T);
    if ~BIND(d).binds
        fprintf(['[NOTA] En esta densidad la mascara NO MUERDE (ver seccion 3b): el servidor\n' ...
                 '       es siempre el satelite mas alto y nunca baja de %.1f deg. Las filas\n' ...
                 '       siguientes deben salir practicamente IDENTICAS; lo unico que cambia\n' ...
                 '       es el numero de interferentes inter-satelite.\n'], min(BIND(d).elev_min));
    end
    for c = 1:nC
        fprintf('\n--- %s ---\n', nm{c});
        fprintf('%8s %8s %9s %10s %9s %9s %9s %9s %8s %8s %9s %10s %s\n', ...
            'minElev','cob.geo','SINRe_p5','cob.calid','margFlr','Rp5_bor','Rp5_glob', ...
            'R_agg','Jain','SE_e_p5','penaliz','elev_med','veredicto');
        for i = 1:nE
            K = R{d}{i}.K{c};  v = K.viab;  dg = R{d}{i}.diag;
            fprintf('%8d %8.2f %9.3f %10.2f %9.3f %9.2f %9.2f %9.2f %8.3f %8.3f %9.2f %10.1f %s\n', ...
                minElevList(i), dg.covFrac, v.SINR_edge_p5, v.coverage, v.margin_floor_dB, ...
                K.edge.R_p5_Mbps, K.all.R_p5_Mbps, K.all.R_agg_Mbps/1e3, K.all.jain, ...
                v.SE_edge_p5, K.penalty_mean, dg.elev_mean, v.verdict);
        end
    end

    fprintf('\n--- Penalizacion por interferencia (C/N - SINR) [dB], media ---\n');
    fprintf('%8s', 'minElev'); for c = 1:nC, fprintf(' %13s', nm{c}); end; fprintf('\n');
    for i = 1:nE
        fprintf('%8d', minElevList(i));
        for c = 1:nC, fprintf(' %13.2f', R{d}{i}.K{c}.penalty_mean); end
        fprintf('\n');
    end

    % Recorrido del KPI en el eje: cuantifica la degeneracion
    v1 = arrayfun(@(i) R{d}{i}.K{1}.viab.SINR_edge_p5, 1:nE);
    fprintf('\nRecorrido de SINR_edge_p5 (reuse1) en TODO el eje: %.3f dB\n', max(v1)-min(v1));
end

%% ------------------------------------------------------------------------
%  7. TAREA 2: COMPRESION ANGULAR MEDIDA
%  ------------------------------------------------------------------------
fprintf('\n\n=======================================================================\n');
fprintf('  TAREA 2: COMPRESION ANGULAR MEDIDA (el mecanismo)\n');
fprintf('=======================================================================\n\n');
compr = cell(1,nD);  sepEf = cell(1,nD);
for d = 1:nD
    fprintf('--- Densidad %s (T=%d) ---\n', densList(d).name, densList(d).T);
    fprintf('%8s %10s %10s %12s %12s %12s %12s %12s\n', ...
        'minElev','elev_med','elev_p5','theta_med','theta_max','sep_ef[deg]','sep/nom','SINRe_p5 r1');
    compr{d} = zeros(1,nE);  sepEf{d} = zeros(1,nE);
    for i = 1:nE
        dg = R{d}{i}.diag;
        sepEf{d}(i) = cell_sep_eff(dg.elev_mean, BL0.spacing_km, ...
                                   cfg.constellations(1).h, cfg.const.Re);
        compr{d}(i) = sepEf{d}(i) / sepNom;
        fprintf('%8d %10.2f %10.2f %12.4f %12.4f %12.4f %12.3f %12.3f\n', ...
            minElevList(i), dg.elev_mean, dg.elev_p5, dg.theta_mean_deg, dg.theta_max_deg, ...
            sepEf{d}(i), compr{d}(i), R{d}{i}.K{1}.viab.SINR_edge_p5);
    end
    fprintf('\n');
end

% Senal operativa de guarda por punto. El criterio (umbral = circunradio de celda
% s/sqrt(3)) vive en grid_outside_cluster, la MISMA funcion que usa el motor en
% ffr_context: aqui no se reimplementa, para que el diagnostico del experimento y el
% aviso del motor no puedan divergir (ver la cabecera de grid_outside_cluster).
[~, ~, thLim0] = grid_outside_cluster(0, BL0.spacing_deg);
fprintf('Guarda operativa (max theta frente al circunradio s/sqrt(3) = %.4f deg, s = %.4f deg):\n', ...
    thLim0, BL0.spacing_deg);
anySesg = false;  pctMax = 0;
for d = 1:nD
    for i = 1:nE
        [isOut, pctC] = grid_outside_cluster(R{d}{i}.diag.theta_max_deg, BL0.spacing_deg);
        pctMax = max(pctMax, pctC);
        if isOut
            fprintf('  %s minElev=%2d: max theta = %.3f deg = %.1f%% del circunradio -> SESGADO\n', ...
                densList(d).name, minElevList(i), R{d}{i}.diag.theta_max_deg, pctC);
            anySesg = true;
        end
    end
end
if ~anySesg
    fprintf(['  Ningun punto supera el circunradio (peor caso %.1f%%): NINGUN PUNTO SESGADO\n' ...
             '  (guarda analitica %.2f anillos, constante en todo el eje).\n'], pctMax, guard);
end

% Correlacion compresion <-> KPI, en la densidad donde la mascara MUERDE
for d = 1:nD
    sr1 = arrayfun(@(i) R{d}{i}.K{1}.viab.SINR_edge_p5, 1:nE);
    srf = arrayfun(@(i) R{d}{i}.K{4}.viab.SINR_edge_p5, 1:nE);
    if (max(sr1)-min(sr1)) < 1e-6
        fprintf('\n[%s] el KPI no varia en el eje (mascara no muerde): la correlacion no aplica.\n', ...
            densList(d).name);
        continue;
    end
    fprintf('\n[%s] Correlacion sobre los %d puntos:\n', densList(d).name, nE);
    fprintf('  corr(sep_efectiva, SINR_edge_p5 reuse1) = %+.4f\n', corr_(compr{d}, sr1));
    fprintf('  corr(sep_efectiva, SINR_edge_p5 ffr(3)) = %+.4f\n', corr_(compr{d}, srf));
end
fprintf(['-> Correlacion POSITIVA y alta = al comprimirse la reticula (sep/nom baja),\n' ...
         '   la SINR de borde baja con ella: la compresion EXPLICA la caida.\n']);

%% ------------------------------------------------------------------------
%  7b. MECANISMO AISLADO: KPIs ESTRATIFICADOS por elevacion REAL del servidor
%  ------------------------------------------------------------------------
%  El barrido de mascara mezcla dos efectos (cambia la cobertura Y la geometria).
%  Aqui se aisla la elevacion: UNA sola simulacion (T=66, mascara 10 deg, el rango de
%  elevacion mas ancho) y los KPIs se calculan por BINS de elevacion del servidor con
%  la mascara de muestras que ya admite compute_kpis. Es el metodo con el que E5
%  midio los ~7 dB, ahora sobre todos los esquemas y con la compresion al lado.
fprintf('\n\n=======================================================================\n');
fprintf('  TAREA 2b: MECANISMO AISLADO (KPIs por BIN de elevacion, T=66, mask 10)\n');
fprintf('=======================================================================\n\n');
STRAT = elevation_strata(cfg, cases, BL0, sepNom);

%% ------------------------------------------------------------------------
%  8. MARGEN DE VIABILIDAD: minElev MINIMA con verdict = 'viable'
%  ------------------------------------------------------------------------
fprintf('\n\n=======================================================================\n');
fprintf('  MARGEN DE VIABILIDAD: minElev MINIMA con verdict = ''viable''\n');
fprintf('=======================================================================\n');
fprintf('(H1 en este eje: hasta que elevacion aguanta cada esquema)\n');

MARGIN = cell(1,nD);
for d = 1:nD
    fprintf('\n--- Densidad %s (T=%d)%s ---\n', densList(d).name, densList(d).T, ...
        repmat('   [MASCARA NO MUERDE: eje degenerado]', 1, ~BIND(d).binds));
    fprintf('%-14s %16s %18s %14s  %s\n','esquema','minElev min viable', ...
        'elev cruce interp','margen 0dB @10','veredictos (45 -> 10 deg)');
    Md = struct('scheme',{},'minElev_min_viable',{},'elev_cross',{},'verdicts',{});
    for c = 1:nC
        ver = cell(1,nE);  mth = nan(1,nE);
        for i = 1:nE
            ver{i} = R{d}{i}.K{c}.viab.verdict;
            mth(i) = R{d}{i}.K{c}.viab.margin_th_dB;
        end
        okv = strcmp(ver,'viable');
        if any(okv), emin = min(minElevList(okv)); else, emin = NaN; end

        % Cruce interpolado de margin_th_dB = 0 (el eje va de mas a menos elevacion)
        ecross = NaN;
        for k = 1:nE-1
            if ~isnan(mth(k)) && ~isnan(mth(k+1)) && mth(k) > 0 && mth(k+1) <= 0
                t = mth(k) / (mth(k) - mth(k+1));
                ecross = minElevList(k) + t*(minElevList(k+1) - minElevList(k));
                break;
            end
        end
        vs = strjoin(cellfun(@(s) s(1:3), ver, 'UniformOutput', false), '>');
        fprintf('%-14s %16s %18s %14.3f  %s\n', nm{c}, num2str(emin,'%d'), ...
            num2str(ecross,'%.1f'), mth(end), vs);

        Md(c).scheme = nm{c};
        Md(c).minElev_min_viable = emin;
        Md(c).elev_cross = ecross;
        Md(c).verdicts = ver;
    end
    MARGIN{d} = Md;

    fprintf('\nGanancia de la FFR sobre reuse1 (SINR_edge_p5, dB):\n');
    fprintf('%8s %12s %12s %12s\n','minElev','reuse1','ffr(3)','ffr(3)-r1');
    for i = 1:nE
        a = R{d}{i}.K{1}.viab.SINR_edge_p5;  b = R{d}{i}.K{4}.viab.SINR_edge_p5;
        fprintf('%8d %12.3f %12.3f %12.3f\n', minElevList(i), a, b, b-a);
    end
end

%% ------------------------------------------------------------------------
%  9. FIGURAS
%  ------------------------------------------------------------------------
if ~exist('figs_e3b','dir'), mkdir('figs_e3b'); end
make_e3b_figures(R, densList, cases, minElevList, elGrid, sepRad, sepTan, sepNom, cfg, STRAT);

%% ------------------------------------------------------------------------
%  10. TAREA 4: contraste de cfg.e3.powerMode (cabo pendiente de E3)
%  ------------------------------------------------------------------------
fprintf('\n\n=======================================================================\n');
fprintf('  TAREA 4: CONTRASTE cfg.e3.powerMode (total_const vs per_beam_const)\n');
fprintf('=======================================================================\n');
fprintf(['PREDICCION: no debe cambiar la conclusion. La potencia entra en Cpsd, que es\n' ...
         'COMUN a la portadora y a la interferencia INTRA (mismo satelite, mismo Cpsd en\n' ...
         'compute_sinr_ffr), luego la C/I es INVARIANTE a la potencia; solo cambia el C/N.\n' ...
         'Como el sistema esta limitado por interferencia (penalizacion C/N-SINR ~ 20-28 dB\n' ...
         '=> I >> N), la SINR apenas debe moverse.\n\n']);

PWR = powermode_contrast(cfg, cases, nWorkers);

%% ------------------------------------------------------------------------
%  11. Guardado
%  ------------------------------------------------------------------------
save('minelev_results.mat', 'res', 'R', 'ptMeta', 'INFO', 'cases', 'nm', ...
     'minElevList', 'densList', 'BIND', 'compr', 'sepEf', 'sepNom', 'elGrid', ...
     'sepRad', 'sepTan', 'MARGIN', 'STRAT', 'PWR', 'MEM', 'cfg', 'alphaStar', ...
     'guard', '-v7.3');
fprintf('\nResultados en minelev_results.mat | figuras en figs_e3b/\n');

%% ========================================================================
%  FUNCIONES LOCALES
%  ========================================================================
function [sepRad_deg, sepTan_deg] = cell_sep_eff(el_deg, s_km, h_km, Re_km)
%CELL_SEP_EFF  Separacion angular efectiva entre centros de celda vecinos.
%   [rad, tan] = cell_sep_eff(el_deg, s_km, h_km)
%
%   Angulo que subtiende EN EL SATELITE un segmento de suelo de longitud s_km,
%   orientado (a) RADIALMENTE, en el plano vertical que contiene al satelite -- la
%   direccion que se escorza -- y (b) TRANSVERSALMENTE (perpendicular). Es la
%   version cuantitativa de la "compresion angular": a nadir vale ~s_km/h, y al
%   bajar la elevacion la radial se comprime aproximadamente como sin(elev).
%
%   Geometria en ENU local con azimut del satelite = Norte (sin perdida de
%   generalidad): el satelite esta en r*[0, cos(el), sin(el)] con r el rango
%   oblicuo, obtenido con Tierra esferica.
%
%   Re_km (OPCIONAL): radio terrestre. CORRECCION (auditoria, hallazgo B10): antes
%   estaba HARDCODEADO a 6371 km (radio medio) mientras TODO el resto del simulador
%   usa cfg.const.Re = 6378.137 km (ecuatorial WGS-84). Ahora los llamadores pasan
%   cfg.const.Re y la constante deja de estar duplicada. Efecto numerico medido:
%   cambio relativo 1e-5 a 3e-4 en la separacion efectiva, que NO altera ninguna de
%   las cifras publicadas de E3b a los 3 decimales con que se reportan
%   (rad/nom 0.767 / 0.525 / 0.209 / 0.053 se mantienen). No toca ningun KPI:
%   cell_sep_eff es una funcion de ANALISIS, no forma parte del pipeline.
if nargin < 4 || isempty(Re_km), Re_km = 6378.137; end
Re = Re_km;
% Rango oblicuo para elevacion el y altitud h (Tierra esferica)
psi  = asind( Re*cosd(el_deg) / (Re + h_km) );        % angulo en el satelite
gam  = 90 - el_deg - psi;                             % angulo central
r    = (Re + h_km) * sind(gam) / cosd(el_deg);        % ley de senos

S  = r * [0, cosd(el_deg), sind(el_deg)];             % satelite en ENU (E,N,U)
A  = [0, 0, 0];                                       % celda de referencia
Br = [0, s_km, 0];                                    % vecina RADIAL (hacia el sat)
Bt = [s_km, 0, 0];                                    % vecina TRANSVERSAL

sepRad_deg = ang_between(A - S, Br - S);
sepTan_deg = ang_between(A - S, Bt - S);
end

% -------------------------------------------------------------------------
function a = ang_between(u, v)
%ANG_BETWEEN  Angulo entre dos vectores [deg].
c = dot(u,v) / (norm(u)*norm(v));
a = acosd(max(-1, min(1, c)));
end

% -------------------------------------------------------------------------
function r = corr_(x, y)
%CORR_  Coeficiente de correlacion de Pearson (sin Statistics Toolbox).
x = x(:) - mean(x(:));  y = y(:) - mean(y(:));
r = (x.'*y) / sqrt((x.'*x)*(y.'*y));
end

% -------------------------------------------------------------------------
function S = elevation_strata(cfgBase, cases, BL0, sepNom)
%ELEVATION_STRATA  KPIs por BIN de elevacion del servidor (mecanismo aislado).
%   Una sola simulacion (T=66, mascara 10 deg) y los KPIs se calculan restringiendo
%   la MUESTRA a cada bin de elevacion, con el argumento `mask` que ya admite
%   compute_kpis. Asi la comparacion entre bins es a poblacion, instantes y
%   clasificacion centro/borde IDENTICOS: lo unico que cambia es la geometria.
cfg = cfgBase;
cfg.constellations(1).T = 66;
cfg.constellations(1).P = 6;
cfg.geom.minElev = 10;                   % rango de elevacion lo mas ancho posible

edges = [10 20 30 40 50 60 90];
nB    = numel(edges) - 1;
nC    = numel(cases);
nm    = cellfun(@(c) c.name, cases, 'UniformOutput', false);

tvec  = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
Nt    = numel(tvec);
sats  = build_constellation(cfg);
users = build_user_grid(cfg);
BL    = build_beam_layout(cfg);
M     = users.M;

nBlk = 20;  if isfield(cfg,'compute') && ~isempty(cfg.compute.timeBlock), nBlk = cfg.compute.timeBlock; end
nBlk = min(nBlk, Nt);  nBlocks = ceil(Nt/nBlk);

acc = repmat(struct('SINR',nan(M,Nt),'CN',nan(M,Nt),'pen',nan(M,Nt), ...
                    'Bu',nan(M,Nt),'isC',false(M,Nt),'alpha',nan(1,Nt)), 1, nC);
elevAcc = nan(M,Nt);  thetaAcc = nan(M,Nt);  covAcc = false(M,Nt);
metaD = cell(1,nC);

for b = 1:nBlocks
    kIdx = ((b-1)*nBlk+1) : min(b*nBlk, Nt);
    tv = tvec(kIdx);
    R_eci = propagate(cfg, sats, tv);  R_ecef = eci2ecef(cfg, R_eci, tv);
    G = compute_geometry(cfg, users, R_ecef);
    Sv = associate_serving(cfg, G);
    L = compute_link_budget(cfg, Sv);
    INTF = compute_interference(cfg, sats, users, G, Sv, L);
    ctx = ffr_context(cfg, sats, users, BL, R_ecef, G, Sv, L, INTF, false);
    clear G R_eci R_ecef INTF Sv L;
    for c = 1:nC
        cc = cfg;
        cc.ffr.scheme = cases{c}.scheme;  cc.ffr.Delta = cases{c}.Delta;
        if ~isnan(cases{c}.alpha), cc.ffr.alpha = cases{c}.alpha; end
        cc.ffr.adaptive = cases{c}.adaptive;
        A = ffr_allocate(cc, ctx, false);
        F = compute_sinr_ffr(cc, ctx, A);
        acc(c).SINR(:,kIdx)=F.SINR_dB; acc(c).CN(:,kIdx)=F.CN_dB; acc(c).pen(:,kIdx)=F.penalty_dB;
        acc(c).Bu(:,kIdx)=A.B_user_Hz;  acc(c).isC(:,kIdx)=A.isCenter; acc(c).alpha(kIdx)=A.alpha;
        if isempty(ctx.SINR_ref_dB), ctx.SINR_ref_dB = A.SINR_ref_dB; end
        if b==1, metaD{c} = struct('Delta',A.Delta); end
    end
    elevAcc(:,kIdx)=ctx.elev_deg; thetaAcc(:,kIdx)=ctx.theta_deg; covAcc(:,kIdx)=ctx.cov;
    clear ctx;
end

fprintf('%10s %8s %10s %12s %12s', 'bin elev','n','theta_med','sep_ef[deg]','sep/nom');
for c = 1:nC, fprintf(' %13s', nm{c}); end
fprintf('\n');

S = struct('edges',edges,'bins',[],'n',[],'theta',[],'sep',[],'compr',[],'SINR',[],'penalty',[],'verdict',{{}});
S.SINR = nan(nB,nC);  S.penalty = nan(nB,nC);  S.verdict = cell(nB,nC);
S.n = zeros(1,nB); S.theta = nan(1,nB); S.sep = nan(1,nB); S.compr = nan(1,nB);
S.bins = zeros(1,nB);
for b = 1:nB
    msk = covAcc & elevAcc >= edges(b) & elevAcc < edges(b+1);
    S.bins(b) = (edges(b)+edges(b+1))/2;
    S.n(b) = sum(msk(:));
    if S.n(b) == 0, fprintf('%6d-%-3d %8d   (sin muestras)\n', edges(b), edges(b+1), 0); continue; end
    S.theta(b) = mean(thetaAcc(msk));
    S.sep(b)   = cell_sep_eff(mean(elevAcc(msk)), BL0.spacing_km, ...
                              cfg.constellations(1).h, cfg.const.Re);
    S.compr(b) = S.sep(b) / sepNom;
    fprintf('%6d-%-3d %8d %10.4f %12.4f %12.3f', edges(b), edges(b+1), S.n(b), ...
        S.theta(b), S.sep(b), S.compr(b));
    for c = 1:nC
        cc = cfg; cc.ffr.scheme = cases{c}.scheme; cc.ffr.Delta = cases{c}.Delta;
        FF = struct('SINR_dB',acc(c).SINR,'CN_dB',acc(c).CN,'penalty_dB',acc(c).pen);
        AL = struct('B_user_Hz',acc(c).Bu,'isCenter',acc(c).isC,'scheme',cc.ffr.scheme, ...
                    'Delta',metaD{c}.Delta,'alpha',acc(c).alpha);
        K = compute_kpis(cc, FF, AL, msk);
        S.SINR(b,c) = K.viab.SINR_edge_p5;  S.penalty(b,c) = K.penalty_mean;
        S.verdict{b,c} = K.viab.verdict;
        fprintf(' %13.3f', S.SINR(b,c));
    end
    fprintf('\n');
end

ok = ~isnan(S.compr) & ~isnan(S.SINR(:,1)).';
fprintf(['\ncorr(compresion, SINR_edge_p5 reuse1) sobre los bins = %+.4f\n' ...
         'Recorrido de SINR_edge_p5 reuse1 entre bins = %.2f dB  <-- el efecto de la elevacion\n'], ...
    corr_(S.compr(ok), S.SINR(ok,1).'), max(S.SINR(ok,1))-min(S.SINR(ok,1)));
end

% -------------------------------------------------------------------------
function make_e3b_figures(R, densList, cases, minElevList, elGrid, sepRad, sepTan, sepNom, cfg, STRAT)
%make_e3b_figures  las cuatro figuras de E3b.
nE = numel(minElevList);  nC = numel(cases);  nD = numel(densList);
col = lines(nC);
nm  = cellfun(@(c) c.name, cases, 'UniformOutput', false);
gth = cfg.viab.gamma_th_dB;  gfl = cfg.viab.gamma_floor_dB;

% (a) SINR_edge_p5 vs minElev, una columna por densidad
figure('Name','E3b(a) SINR de borde vs minElev','Position',[60 60 1150 480]);
for d = 1:nD
    subplot(1,nD,d); hold on; grid on;
    for c = 1:nC
        plot(minElevList, arrayfun(@(i) R{d}{i}.K{c}.viab.SINR_edge_p5, 1:nE), '-o', ...
            'Color', col(c,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(c,:), 'MarkerSize', 4);
    end
    yline(gth, 'k--', sprintf('\\gamma_{th} = %.1f', gth), 'LineWidth', 1.2);
    yline(gfl, 'r--', sprintf('\\gamma_{floor} = %.1f', gfl), 'LineWidth', 1.2);
    set(gca,'XDir','reverse');
    xlabel('minElev [deg]'); ylabel('SINR_{edge,p5} [dB]');
    title(sprintf('%s (T=%d)', densList(d).name, densList(d).T));
    if d == 1, legend(nm,'Location','southwest','FontSize',7); end
end
sgtitle('E3b(a): SINR de borde vs mascara de elevacion');
save_fig(gcf, 'figs_e3b', 'e3b_a_sinr_edge_vs_minelev');

% (b) Mapa de veredictos
figure('Name','E3b(b) Veredicto vs minElev','Position',[60 60 1150 460]);
lev = containers.Map({'inviable','marginal','viable','sin-datos'}, {1,2,3,0});
for d = 1:nD
    subplot(1,nD,d);
    Mv = zeros(nC, nE);
    for c = 1:nC, for i = 1:nE, Mv(c,i) = lev(R{d}{i}.K{c}.viab.verdict); end, end
    imagesc(Mv, [0 3]); colormap([0.85 0.85 0.85; 0.80 0.20 0.20; 0.95 0.75 0.20; 0.20 0.65 0.30]);
    set(gca,'YTick',1:nC,'YTickLabel',nm,'FontSize',8);
    set(gca,'XTick',1:nE,'XTickLabel',arrayfun(@(x) sprintf('%d',x), minElevList, 'UniformOutput', false));
    xlabel('minElev [deg]'); title(sprintf('%s (T=%d)', densList(d).name, densList(d).T));
    for c = 1:nC
        for i = 1:nE
            text(i, c, R{d}{i}.K{c}.viab.verdict(1:3), 'HorizontalAlignment','center', ...
                'FontSize',6,'Color','w','FontWeight','bold');
        end
    end
end
sgtitle('E3b(b): VEREDICTO (rojo=inviable, ambar=marginal, verde=viable)');
save_fig(gcf, 'figs_e3b', 'e3b_b_verdict_map');

% (c) Penalizacion vs minElev
figure('Name','E3b(c) Penalizacion vs minElev','Position',[60 60 1150 480]);
for d = 1:nD
    subplot(1,nD,d); hold on; grid on;
    for c = 1:nC
        plot(minElevList, arrayfun(@(i) R{d}{i}.K{c}.penalty_mean, 1:nE), '-o', ...
            'Color', col(c,:), 'LineWidth', 1.6, 'MarkerSize', 4);
    end
    set(gca,'XDir','reverse');
    xlabel('minElev [deg]'); ylabel('C/N - SINR [dB]');
    title(sprintf('%s (T=%d)', densList(d).name, densList(d).T));
    if d == 1, legend(nm,'Location','best','FontSize',7); end
end
sgtitle('E3b(c): penalizacion por interferencia vs minElev');
save_fig(gcf, 'figs_e3b', 'e3b_c_penalty');

% (d) EL MECANISMO: compresion analitica + KPI estratificado por elevacion REAL
figure('Name','E3b(d) Compresion angular','Position',[60 60 1150 460]);
subplot(1,2,1); hold on; grid on;
plot(elGrid, sepRad, 'LineWidth', 1.8);
plot(elGrid, sepTan, '--', 'LineWidth', 1.4);
yline(sepNom, ':', 'nominal s', 'LineWidth', 1.2);
xlabel('elevacion del servidor [deg]'); ylabel('separacion efectiva [deg]');
title('Separacion angular entre celdas vecinas');
legend({'radial (se escorza)','transversal','nominal'},'Location','southeast','FontSize',8);

subplot(1,2,2);
ok = ~isnan(STRAT.compr);
yyaxis left; hold on; grid on;
plot(STRAT.bins(ok), STRAT.compr(ok), '-o', 'LineWidth', 1.8, 'MarkerSize', 5);
ylabel('sep. efectiva / nominal');
yyaxis right;
plot(STRAT.bins(ok), STRAT.SINR(ok,1), '-s', 'LineWidth', 1.8, 'MarkerSize', 5);
plot(STRAT.bins(ok), STRAT.SINR(ok,4), '-^', 'LineWidth', 1.8, 'MarkerSize', 5);
ylabel('SINR_{edge,p5} [dB]');
xlabel('elevacion REAL del servidor [deg]');
title('Compresion y KPI, por bin de elevacion');
legend({'compresion','SINR reuse1','SINR ffr(3)'},'Location','northwest','FontSize',8);
sgtitle('E3b(d): compresion angular = el MECANISMO (KPI estratificado, no barrido de mascara)');
save_fig(gcf, 'figs_e3b', 'e3b_d_compression');

fprintf('\n4 figuras exportadas a figs_e3b/\n');
end

% -------------------------------------------------------------------------
function PWR = powermode_contrast(cfgBase, cases, nWorkers)
%POWERMODE_CONTRAST  Repite 3 puntos de E3 con las dos hipotesis de potencia.
%   Mismo acoplamiento que run_e3_beamdensity. Al estrechar el haz la directividad
%   sube 20log10(th_ref/th); lo que decide el modo es si la potencia POR HAZ baja en
%   el mismo factor (potencia total del satelite constante) o no:
%     total_const     -> EIRPdens_eff = EIRPdens_ref                  (se cancelan)
%     per_beam_const  -> EIRPdens_eff = EIRPdens_ref + 20log10(th_ref/th)
%
%   EL ACOPLAMIENTO VIAJA POR EIRPdensity_dBWMHz, NO POR Gmax_dBi (auditoria, G2).
%   cfg.radio.Gmax_dBi NO entra en la fisica: la densidad EIRP ya INCLUYE la ganancia
%   de pico del haz por definicion de EIRP, y sumarla aparte seria doble
%   contabilizacion. Antes se sobrescribia aqui Gmax_dBi = Gmax_ref + dG, lo cual era
%   un NO-OP (no cambiaba ningun numero) pero sugeria un acoplamiento que el codigo no
%   implementa. Se ha eliminado. El aviso esta en la cabecera de compute_link_budget.
cfg = cfgBase;
cfg.geom.minElev = 25;                       % el valor nominal del perfil de estudio
thetas   = [2.2 1.5 0.7];                    % 3 puntos de E3
th_ref   = 1.5;
EIRPref  = cfg.radio.EIRPdensity_dBWMHz;
modes    = {'total_const','per_beam_const'};

% Guarda preservada (variante B de E3): nRings para que el cluster cubra la rejilla
nRb = zeros(size(thetas));
for i = 1:numel(thetas)
    sdeg = rad2deg(sqrt(3)*sin(deg2rad(thetas(i))/2));
    skm  = cfg.constellations(1).h * tand(sdeg);
    n = 2; while (n*skm*sqrt(3)/2 - cfg.ground.radius_km)/skm < 1 && n < 40, n = n + 1; end
    nRb(i) = n;
end

pts = {};  meta = [];
for m = 1:numel(modes)
    for i = 1:numel(thetas)
        th = thetas(i);
        dG = 20*log10(th_ref/th);
        if strcmp(modes{m},'total_const'), eirp = EIRPref; else, eirp = EIRPref + dG; end
        pts{end+1} = struct( ...
            'cfg_over', struct('radio', struct('beamwidth3dB_deg', th, ...
                                               'EIRPdensity_dBWMHz', eirp), ...
                               'beams', struct('nRings', nRb(i), 'nBeams', [])), ...
            'cases', {cases}, 'label', sprintf('%s-th%.1f', modes{m}, th)); %#ok<AGROW>
        meta(end+1,:) = [m, i]; %#ok<AGROW>
    end
end

fprintf('Ejecutando %d puntos (%d thetas x %d modos)...\n', numel(pts), numel(thetas), numel(modes));
r = run_sweep_points(cfg, pts, nWorkers, struct('verbose',false));

fprintf('\n%8s %14s %14s %12s %12s %12s %12s\n', ...
    'HPBW','SINRe total','SINRe perbeam','dif[dB]','C/N tot','C/N pb','penaliz tot');
PWR = struct('theta',{},'SINR_total',{},'SINR_perbeam',{},'dSINR',{},'verdict_total',{},'verdict_perbeam',{});
maxd = 0;
for i = 1:numel(thetas)
    kT = find(meta(:,1)==1 & meta(:,2)==i);
    kP = find(meta(:,1)==2 & meta(:,2)==i);
    vT = r{kT}.K{4}.viab;  vP = r{kP}.K{4}.viab;    % ffr(3,a*) como referencia
    cnT = r{kT}.K{4}.penalty_mean;  cnP = r{kP}.K{4}.penalty_mean;
    d   = vP.SINR_edge_p5 - vT.SINR_edge_p5;
    maxd = max(maxd, abs(d));
    fprintf('%8.2f %14.3f %14.3f %12.4f %12.2f %12.2f %12.2f\n', ...
        thetas(i), vT.SINR_edge_p5, vP.SINR_edge_p5, d, cnT, cnP, cnT);
    PWR(i).theta = thetas(i);
    PWR(i).SINR_total = vT.SINR_edge_p5;
    PWR(i).SINR_perbeam = vP.SINR_edge_p5;
    PWR(i).dSINR = d;
    PWR(i).verdict_total = vT.verdict;
    PWR(i).verdict_perbeam = vP.verdict;
end

sameVerdict = all(arrayfun(@(k) strcmp(PWR(k).verdict_total, PWR(k).verdict_perbeam), 1:numel(PWR)));
fprintf('\nMaxima diferencia de SINR_edge_p5 entre modos: %.4f dB\n', maxd);
fprintf('Veredictos identicos en los %d puntos: %s\n', numel(PWR), string(sameVerdict));
if maxd < 0.5 && sameVerdict
    fprintf(['VEREDICTO: la conclusion de E3 NO cambia con powerMode. Confirmada la\n' ...
             'prediccion: la potencia es comun a C y a I intra, luego la C/I es invariante,\n' ...
             'y el sistema esta limitado por interferencia (I >> N).\n']);
else
    fprintf(['[OBJECION] La conclusion SI depende de powerMode (max %.3f dB, veredictos\n' ...
             'iguales = %s). Hay que revisar E3: el sistema no estaria tan limitado por\n' ...
             'interferencia como se supuso.\n'], maxd, string(sameVerdict));
end
end

