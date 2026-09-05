%% RUN_CONVERGENCE_NRINGS  Estudio de convergencia del tamano del cluster de haces.
%
%  DOBLE PROPOSITO:
%   (a) VALIDAR que cfg.beams.nRings = 5 (91 haces) captura toda la interferencia
%       INTRA-satelite relevante, es decir que anadir mas anillos ya no mueve el KPI.
%   (b) DECIDIR si nRings sirve como EJE DE SATURACION secundario en E3 o es solo
%       precision del modelo. La FASE A dejo la pregunta abierta: la densidad ORBITAL
%       (T) no satura el sistema porque la interferencia INTER-satelite es
%       despreciable (I_inter/N <= -21.7 dB con 4000 satelites), luego el unico eje
%       que escala la interferencia DOMINANTE (la intra) es la densidad de HACES
%       CO-CANAL sobre la misma zona. Este script mide si ese eje tiene recorrido.
%
%  PREGUNTA FISICA: el patron de haz (Bessel 38.811 + suelo de lobulos) atenua los
%  haces lejanos, asi que la potencia interferente que aporta el anillo n decae con
%  n. Si a partir de cierto n la aportacion marginal es despreciable, anadir anillos
%  es PRECISION DEL MODELO y no saturacion. Referencia externa: en redes celulares y
%  VLC la practica estandar son layouts de 3 anillos ("three-tier"), cf. Genoves
%  Guzman et al., IEEE Trans. Commun. 68(10), 2020.
%
%  METODOLOGIA (importante para la defensa: aqui NO hay fisica nueva)
%  ------------------------------------------------------------------
%  El pipeline es EXACTAMENTE el de run_one_density (mismo orden de llamadas al
%  motor, mismo troceado temporal). Lo unico que anade este script es una
%  DESCOMPOSICION POR ANILLO de la interferencia intra-satelite: se reparte la MISMA
%  suma co-canal que calcula compute_sinr_ffr
%
%      I_intra_lin = sum_{b in co-canal} Glin(m,b,k)  -  Glin(m,b_propio,k)
%
%  entre los anillos hexagonales de BL.ring, respetando la regla de co-canal del
%  motor (subband = 0 -> todos los haces; subband = c -> solo los de color c). Que
%  es una particion EXACTA y no un modelo aparte se COMPRUEBA: la suma de las
%  aportaciones por anillo se contrasta contra el I_intra_dBW que devuelve
%  compute_sinr_ffr y se aborta si difieren.
%
%  GUARDA (sesgo que hay que declarar): la rejilla de usuarios tiene radio 40 km y
%  las celdas separacion s ~ 12.5 km, luego con nRings bajo el cluster NO cubre la
%  rejilla: los usuarios del borde caen fuera del ultimo anillo, se asocian a un haz
%  muy desapuntado y ademas ven MENOS interferentes de los que les tocarian. En esos
%  casos la comparacion esta SESGADA y el script lo dice explicitamente.
%
%  Uso: ejecutar en la raiz del proyecto. Resultados en convergence_nrings.mat.
%  NO modifica el motor: E0 y run_ffr_demo quedan bit-identicos (este fichero solo
%  se anade).

clear; clc; close all;

fprintf('=======================================================================\n');
fprintf('  ESTUDIO DE CONVERGENCIA DE nRings (tamano del cluster de haces)\n');
fprintf('=======================================================================\n\n');

%% ------------------------------------------------------------------------
%  1. CONFIGURACION: la CONVERGIDA de la FASE A, con todo lo demas FIJO
%  ------------------------------------------------------------------------
cfg = config_default();

% Densidad orbital del punto de referencia de la FASE A
cfg.constellations(1).T = 1584;
cfg.constellations(1).P = 72;

% Rejilla y ventana de la FASE A: step_km = 4 -> 317 usuarios, Nt = 61.
% PRECISION (no confundir dos cosas distintas):
%   - Nt = 61 SI es el valor CONVERGIDO por el criterio (tol 0.2 dB).
%   - step_km = 4 NO es el convergido: el criterio converge en 6 km
%     (step_conv = 6 en convergence_studyA.mat). El 4 es el valor DE TRABAJO,
%     mas FINO que el convergido, adoptado por CONSERVADURISMO y por coherencia
%     con los experimentos ya ejecutados (E3, E3b, E4, EB y este estudio), no
%     por convergencia.
% Solo es un comentario: no interviene en ningun calculo.
cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 4;
cfg.time.dt          = 120;
cfg.time.duration    = 7200;

% Umbral centro/borde por CUANTIL (mismo criterio que run_ffr_demo / E5)
cfg.ffr.tau_mode = 'quantile';
cfg.ffr.tau_q    = 50;

% Troceado temporal: nBeams alto sube el termino M*nBeams*Ntb del modelo de
% memoria, y ademas la maquina de desarrollo tiene poca RAM libre. Con bloques de
% 10 instantes el pico se queda por debajo de 0.4 GB incluso con 169 haces.
% El troceado es EXACTO (test_timeblock_invariance: max|dif| = 0).
timeBlock = 10;

% Barrido del cluster. TODO LO DEMAS SE MANTIENE FIJO: spacing y HPBW no se tocan,
% asi que anadir anillos SOLO anade haces co-canal mas lejanos.
nRingsList = [2 3 4 5 6 7];
nBeamsOf   = @(n) 3*n.*(n+1) + 1;          % 19, 37, 61, 91, 127, 169

% Esquemas evaluados (comparten ctx, poblacion e instantes)
cases = { ...
    struct('scheme','reuse1','Delta',1,'alpha',NaN,'name','reuse1'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',0.4,'name','ffr(D=3,a=0.4)') };
nC = numel(cases);

fprintf('Perfil (FASE A convergida): T=%d P=%d | rejilla %.0f km paso %.0f km | dt=%.0f s\n', ...
    cfg.constellations(1).T, cfg.constellations(1).P, cfg.ground.radius_km, ...
    cfg.ground.step_km, cfg.time.dt);
fprintf('HPBW=%.2f deg y separacion inter-haz SIN CAMBIAR | tau por cuantil %.0f%% | timeBlock=%d\n', ...
    cfg.radio.beamwidth3dB_deg, cfg.ffr.tau_q, timeBlock);
fprintf('Barrido nRings = [%s]  ->  nBeams = [%s]\n\n', ...
    num2str(nRingsList), num2str(nBeamsOf(nRingsList)));

%% ------------------------------------------------------------------------
%  2. Piezas independientes del cluster (se construyen UNA vez)
%  ------------------------------------------------------------------------
tvec  = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
Nt    = numel(tvec);
sats  = build_constellation(cfg);
users = build_user_grid(cfg);
M     = users.M;
N     = sats.N;

fprintf('Escenario: M=%d usuarios | Nt=%d instantes | N=%d satelites\n\n', M, Nt, N);

nP  = numel(nRingsList);
RES = cell(1, nP);

% PRECALENTAMIENTO de la cache P.618 (atm_loss_dB) FUERA del cronometro. Sin esto
% el PRIMER punto del barrido carga con la construccion de la tabla (~80 s) y la
% tabla de coste de la TAREA 3 seria ilegible: el punto 1 pareceria ~300x mas caro
% que los demas por un coste fijo que NO depende de nRings. Ninguno de los 5 campos
% de la clave de la cache depende del cluster, luego una sola tabla sirve para todo
% el barrido. No afecta a ningun numero (test_atm_cache: identico bit a bit).
fprintf('[warmup] construyendo la tabla ITU-R P.618 fuera del cronometro...\n');
tW = tic;
atm_loss_dB(cfg, 45);
tWarm = toc(tW);
fprintf('[warmup] listo (%.1f s)\n\n', tWarm);

% CRONOMETRO DE EXTREMO A EXTREMO del estudio, arrancado
% DESPUES del warmup por el mismo motivo que la tabla de coste de la TAREA 3.
tStudy = tic;

% Baseline de memoria del proceso, medido DESPUES del warmup: asi el pico
% muestreado de cada punto es el INCREMENTO atribuible al punto y no arrastra la
% tabla P.618 ni la huella base de MATLAB.
memBase = 0;
try
    mm = memory; memBase = mm.MemUsedMATLAB / 2^30;
catch
end

%% ------------------------------------------------------------------------
%  3. BARRIDO DE nRings
%  ------------------------------------------------------------------------
for p = 1:nP
    nR = nRingsList(p);
    cfgP = cfg;
    cfgP.beams.nRings = nR;
    cfgP.beams.nBeams = [];                 % nRings manda

    BL = build_beam_layout(cfgP);
    nB = BL.nBeams;

    fprintf('-----------------------------------------------------------------------\n');
    fprintf('nRings = %d  (%d haces)\n', nR, nB);

    % --- GUARDA: cuanto sobresale la rejilla del cluster ---------------------
    dCell   = sqrt(sum(BL.cell_km.^2, 2));            % distancia de cada celda al centro
    dOutMin = min(dCell(BL.ring == nR));              % anillo externo, caso conservador
    dOutMax = max(dCell(BL.ring == nR));
    guardMin = (dOutMin - cfgP.ground.radius_km) / BL.spacing_km;   % conservador
    guardMax = (dOutMax - cfgP.ground.radius_km) / BL.spacing_km;

    fprintf('  s = %.4f deg = %.2f km | anillo externo a %.1f-%.1f km | rejilla %.0f km\n', ...
        BL.spacing_deg, BL.spacing_km, dOutMin, dOutMax, cfgP.ground.radius_km);
    if guardMin >= 1
        fprintf('  GUARDA: %.2f anillos (conservador) / %.2f (esquinas)  -> suficiente\n', ...
            guardMin, guardMax);
    elseif guardMin > 0
        fprintf('  GUARDA: %.2f anillos (conservador) / %.2f (esquinas)  -> ESCASA (<1 anillo)\n', ...
            guardMin, guardMax);
    else
        fprintf(['  GUARDA: %.2f anillos -> SIN GUARDA: la rejilla SOBRESALE del cluster.\n' ...
                 '          Los usuarios del borde se asocian a un haz muy desapuntado y ven\n' ...
                 '          MENOS interferentes de los que les corresponden: PUNTO SESGADO.\n'], guardMin);
    end

    %% --- Bucle de BLOQUES TEMPORALES (mismo pipeline que run_one_density) ---
    nBlk    = min(timeBlock, Nt);
    nBlocks = ceil(Nt / nBlk);

    % Acumuladores de MUESTRAS [M x Nt] (no de KPIs parciales: los percentiles se
    % calculan al final sobre la matriz completa, igual que en run_one_density).
    acc = repmat(struct('SINR',nan(M,Nt),'CN',nan(M,Nt),'pen',nan(M,Nt), ...
                        'Bu',nan(M,Nt),'isC',false(M,Nt),'alpha',nan(1,Nt)), 1, nC);
    covAcc  = false(M,Nt);
    elevAcc = nan(M,Nt);
    thetaAcc= nan(M,Nt);

    % Acumuladores de la DESCOMPOSICION POR ANILLO: potencia interferente lineal
    % agregada (suma sobre muestras cubiertas) de cada anillo 0..nR, por esquema.
    Pring   = zeros(nC, nR+1);
    Ptot    = zeros(1, nC);      % control: potencia intra total segun el motor
    nCovS   = zeros(1, nC);
    decErr  = 0;                 % max error relativo de la descomposicion

    metaC   = cell(1,nC);
    memPeakSampled = memBase;
    tPoint  = tic;

    for b = 1:nBlocks
        kIdx = ((b-1)*nBlk + 1) : min(b*nBlk, Nt);
        tv   = tvec(kIdx);

        R_eci  = propagate(cfgP, sats, tv);
        R_ecef = eci2ecef(cfgP, R_eci, tv);
        G      = compute_geometry(cfgP, users, R_ecef);
        S      = associate_serving(cfgP, G);
        L      = compute_link_budget(cfgP, S);
        INTF   = compute_interference(cfgP, sats, users, G, S, L);
        ctx    = ffr_context(cfgP, sats, users, BL, R_ecef, G, S, L, INTF, false);

        % Liberar aqui es lo que ACOTA el pico (si G siguiera viva, el bloque
        % siguiente sumaria su geometria a la anterior).
        clear G R_eci R_ecef INTF S L;

        try
            mm = memory;
            memPeakSampled = max(memPeakSampled, mm.MemUsedMATLAB / 2^30);
        catch
        end

        ctxB = ctx;                      % ctx del bloque (SINR_ref se cachea dentro)
        for c = 1:nC
            cfgC = cfgP;
            cfgC.ffr.scheme = cases{c}.scheme;
            cfgC.ffr.Delta  = cases{c}.Delta;
            if ~isnan(cases{c}.alpha), cfgC.ffr.alpha = cases{c}.alpha; end

            A = ffr_allocate(cfgC, ctxB, false);
            F = compute_sinr_ffr(cfgC, ctxB, A);

            acc(c).SINR(:,kIdx) = F.SINR_dB;
            acc(c).CN(:,kIdx)   = F.CN_dB;
            acc(c).pen(:,kIdx)  = F.penalty_dB;
            acc(c).Bu(:,kIdx)   = A.B_user_Hz;
            acc(c).isC(:,kIdx)  = A.isCenter;
            acc(c).alpha(kIdx)  = A.alpha;

            if isempty(ctxB.SINR_ref_dB), ctxB.SINR_ref_dB = A.SINR_ref_dB; end
            if b == 1, metaC{c} = struct('Delta', A.Delta, 'sched', A.sched); end

            % --- DESCOMPOSICION POR ANILLO (misma suma, repartida) -----------
            [Pr, Pt, nCv, er] = ring_decomposition(ctxB, A, F, BL);
            Pring(c,:) = Pring(c,:) + Pr;
            Ptot(c)    = Ptot(c)    + Pt;
            nCovS(c)   = nCovS(c)   + nCv;
            decErr     = max(decErr, er);
        end

        covAcc(:,kIdx)   = ctxB.cov;
        elevAcc(:,kIdx)  = ctxB.elev_deg;
        thetaAcc(:,kIdx) = ctxB.theta_deg;
        clear ctx ctxB;
    end

    elapsed = toc(tPoint);

    %% --- KPIs sobre las muestras COMPLETAS ------------------------------
    K = cell(1,nC);
    isCenter_ref = [];
    for c = 1:nC
        cfgC = cfgP;
        cfgC.ffr.scheme = cases{c}.scheme;
        cfgC.ffr.Delta  = cases{c}.Delta;
        FF = struct('SINR_dB', acc(c).SINR, 'CN_dB', acc(c).CN, 'penalty_dB', acc(c).pen);
        AL = struct('B_user_Hz', acc(c).Bu, 'isCenter', acc(c).isC, ...
                    'scheme', cfgC.ffr.scheme, 'Delta', metaC{c}.Delta, ...
                    'alpha', acc(c).alpha);
        K{c} = compute_kpis(cfgC, FF, AL);

        % Guarda: la poblacion de borde debe ser la misma en los dos esquemas.
        if isempty(isCenter_ref)
            isCenter_ref = acc(c).isC;
        elseif ~isequal(acc(c).isC, isCenter_ref)
            error('run_convergence_nrings:edgePop', ...
                'La clasificacion centro/borde difiere entre esquemas (nRings=%d).', nR);
        end
    end

    % Validacion de la descomposicion (aborta si no es una particion exacta)
    if decErr > 1e-9
        error('run_convergence_nrings:decomp', ...
            ['La descomposicion por anillo NO reproduce el I_intra del motor ' ...
             '(error relativo %.3g): revisar ring_decomposition.'], decErr);
    end

    %% --- Registro del punto --------------------------------------------
    R = struct();
    R.nRings   = nR;
    R.nBeams   = nB;
    R.guardMin = guardMin;
    R.guardMax = guardMax;
    R.biased   = guardMin <= 0;
    R.spacing_km  = BL.spacing_km;
    R.theta_max   = max(thetaAcc(covAcc));
    R.spacing_deg = BL.spacing_deg;
    % Senal OPERATIVA de "rejilla fuera del cluster". Criterio en la funcion
    % compartida grid_outside_cluster (umbral = circunradio de celda s/sqrt(3)), la
    % misma que usa el motor en ffr_context: aqui no se reimplementa. R.gridOutPct
    % es el offset maximo en % del circunradio, que es lo que separa limpiamente los
    % puntos sanos (92-98%) de los sesgados (118-904%).
    [R.gridOut, R.gridOutPct, R.gridOutLim_deg] = ...
        grid_outside_cluster(R.theta_max, BL.spacing_deg);
    R.elapsed_s   = elapsed;
    R.peak_model_GB = mem_peak_model_GB(M, N, nBlk, nB);
    R.peak_full_GB  = mem_peak_model_GB(M, N, Nt,   nB);
    R.peak_sampled_GB = memPeakSampled - memBase;   % INCREMENTO sobre la linea base
    R.nBlocks   = nBlocks;
    R.covFrac   = sum(covAcc(:)) / numel(covAcc);
    R.elev_mean = mean(elevAcc(covAcc));

    R.scheme = cell(1,nC);
    for c = 1:nC
        Pr   = Pring(c,:) / max(nCovS(c),1);       % potencia media por anillo [W]
        Ptl  = sum(Pr);
        s = struct();
        s.name          = cases{c}.name;
        s.SINR_edge_p5  = K{c}.viab.SINR_edge_p5;
        s.SINR_p5       = K{c}.all.SINR_p5;
        s.SE_edge_p5    = K{c}.viab.SE_edge_p5;
        s.R_p5_edge     = K{c}.edge.R_p5_Mbps;
        s.R_agg_Gbps    = K{c}.all.R_agg_Mbps / 1e3;
        s.verdict       = K{c}.viab.verdict;
        s.I_intra_dBW   = 10*log10(Ptl);                % potencia intra AGREGADA
        s.P_ring        = Pr;                           % [1 x nR+1] media por anillo
        s.pct_ring      = 100 * Pr / Ptl;
        % Aportacion MARGINAL del ultimo anillo (la magnitud clave):
        s.marg_pct      = 100 * Pr(end) / Ptl;
        s.marg_dB       = 10*log10( Ptl / max(Ptl - Pr(end), realmin) );
        R.scheme{c}     = s;
    end
    % Numero de haces co-canal efectivos (regla del motor, no una formula aparte)
    R.nCo_reuse1 = nB - 1;
    R.nCo_ffrEdge = numel(find(ffr_coloring(BL, 3, false) == 1)) - 1;

    RES{p} = R;

    if R.gridOut, outStr = '  -> REJILLA FUERA DEL CLUSTER'; else, outStr = ''; end
    fprintf(['  cobertura %.2f | elev media %.1f deg | max theta %.3f deg ' ...
                '(circunradio %.3f deg -> %.1f%%)%s\n'], ...
        R.covFrac, R.elev_mean, R.theta_max, R.gridOutLim_deg, R.gridOutPct, outStr);
    for c = 1:nC
        s = R.scheme{c};
        fprintf('  %-16s SINR_edge_p5 = %+7.3f dB | SINR_p5 = %+7.3f dB | I_intra = %.2f dBW | ultimo anillo %.2f%% (%.3f dB)\n', ...
            s.name, s.SINR_edge_p5, s.SINR_p5, s.I_intra_dBW, s.marg_pct, s.marg_dB);
    end
    fprintf('  tiempo %.1f s | pico modelo %.2f GB (bloque) / %.2f GB (ventana) | muestreado %.2f GB\n\n', ...
        R.elapsed_s, R.peak_model_GB, R.peak_full_GB, R.peak_sampled_GB);
end

%% ------------------------------------------------------------------------
%  4. TAREA 1 - TABLAS
%  ------------------------------------------------------------------------
fprintf('\n=======================================================================\n');
fprintf('  TAREA 1: BARRIDO DE nRings\n');
fprintf('=======================================================================\n\n');

for c = 1:nC
    fprintf('--- Esquema: %s ---\n', cases{c}.name);
    fprintf('%7s %8s %8s %14s %12s %12s %14s %12s %8s\n', ...
        'nRings','nBeams','guarda','SINR_edge_p5','SINR_p5','I_intra','anillo n [%]','anillo n[dB]','sesgo');
    fprintf('%s\n', repmat('-',1,108));
    for p = 1:nP
        R = RES{p};  s = R.scheme{c};
        if R.biased, flag = 'SESGADO'; else, flag = '-'; end
        fprintf('%7d %8d %8.2f %14.3f %12.3f %12.2f %14.3f %12.4f %8s\n', ...
            R.nRings, R.nBeams, R.guardMin, s.SINR_edge_p5, s.SINR_p5, ...
            s.I_intra_dBW, s.marg_pct, s.marg_dB, flag);
    end
    fprintf('\n');
end

fprintf('Haces co-canal efectivos: reuse1 = nBeams-1 ; ffr(D=3) borde = haces del mismo color - 1\n');
fprintf('%7s %14s %16s\n','nRings','co-canal r1','co-canal borde');
for p = 1:nP
    fprintf('%7d %14d %16d\n', RES{p}.nRings, RES{p}.nCo_reuse1, RES{p}.nCo_ffrEdge);
end

% Perfil COMPLETO por anillo en el cluster mas grande: es la evidencia mas directa
% del decaimiento (una sola geometria, todos los anillos en la misma escala).
pMax = nP;
fprintf('\n--- Perfil de interferencia por ANILLO (nRings = %d, geometria unica) ---\n', RES{pMax}.nRings);
fprintf('%7s', 'anillo'); fprintf(' %12s %12s %12s', 'haces','% del total','acumulado %');
fprintf('\n%s\n', repmat('-',1,52));
sR = RES{pMax}.scheme{1};                 % reuse1: todos los haces son co-canal
cfgMax = cfg;
cfgMax.beams.nRings = RES{pMax}.nRings;
cfgMax.beams.nBeams = [];
BLmax  = build_beam_layout(cfgMax);
cum = 0;
for r = 0:RES{pMax}.nRings
    nb_r = sum(BLmax.ring == r);
    cum  = cum + sR.pct_ring(r+1);
    fprintf('%7d %12d %12.3f %12.3f\n', r, nb_r, sR.pct_ring(r+1), cum);
end

%% ------------------------------------------------------------------------
%  5. TAREA 2 - VEREDICTO
%  ------------------------------------------------------------------------
fprintf('\n\n=======================================================================\n');
fprintf('  TAREA 2: VEREDICTO\n');
fprintf('=======================================================================\n');

tolConv  = 0.2;     % dB                                              % <-- DECISION
tolMarg  = 1.0;     % %  de la potencia interferente total            % <-- DECISION
tolRange = 1.0;     % dB de recorrido del KPI para llamarlo "eje"     % <-- DECISION

% (a) CONVERGENCIA. Criterio de la FASE A: no basta comparar con el SIGUIENTE
%     valor (eso genera falsos positivos si dos puntos coinciden por casualidad);
%     se exige que la diferencia sea < tol contra TODOS los nRings MAYORES.
fprintf('\n(a) CONVERGENCIA  (criterio: |SINR_edge_p5(n) - SINR_edge_p5(n'')| < %.1f dB para TODO n'' > n)\n\n', tolConv);

nmin = nan(1,nC);
for c = 1:nC
    v = arrayfun(@(p) RES{p}.scheme{c}.SINR_edge_p5, 1:nP);
    bi = arrayfun(@(p) RES{p}.biased, 1:nP);

    fprintf('  %s:\n', cases{c}.name);
    fprintf('  %7s %14s %14s %16s\n','nRings','SINR_edge_p5','dif vs n+1','max dif vs n''>n');
    for p = 1:nP
        if p < nP
            d1 = v(p+1) - v(p);
            dm = max(abs(v(p+1:end) - v(p)));
            if bi(p), bStr = '   (punto SESGADO: sin guarda)'; else, bStr = ''; end
            fprintf('  %7d %14.3f %14.3f %16.3f%s\n', nRingsList(p), v(p), d1, dm, bStr);
        else
            fprintf('  %7d %14.3f %14s %16s\n', nRingsList(p), v(p), '-', '-');
        end
    end

    % nRings minimo suficiente: el menor n NO SESGADO que cumple el criterio
    nm = NaN;
    for p = 1:nP-1
        if bi(p), continue; end
        if max(abs(v(p+1:end) - v(p))) < tolConv, nm = nRingsList(p); break; end
    end
    nmin(c) = nm;
    if isnan(nm)
        fprintf('  -> NO converge dentro del rango barrido con tol = %.1f dB.\n\n', tolConv);
    else
        fprintf('  -> nRings MINIMO SUFICIENTE = %d\n\n', nm);
    end
end

% Si ALGUN esquema no converge, no se puede declarar un minimo global (max()
% ignoraria el NaN y daria un falso "converge").
if any(isnan(nmin)), nminGlobal = NaN; else, nminGlobal = max(nmin); end

if isnan(nminGlobal)
    fprintf('  CONCLUSION (a): algun esquema NO converge en el rango barrido.\n');
    fprintf('                  nRings = 5 NO puede declararse suficiente con este barrido.\n');
elseif nminGlobal < 5
    fprintf(['                  nRings = 5 es SUFICIENTE y de hecho EXCESIVO: bastaria con %d\n' ...
             '                  (%d haces en vez de 91). Mantenerlo en 5 es margen de seguridad.\n'], ...
        nminGlobal, nBeamsOf(nminGlobal));
elseif nminGlobal == 5
    fprintf('                  nRings = 5 es EXACTAMENTE el minimo suficiente.\n');
else
    fprintf('                  nRings = 5 se queda CORTO: hace falta %d.\n', nminGlobal);
end

% (b) EJE DE SATURACION. TRES sub-criterios, y hacen falta los tres: un eje de
%     saturacion util tiene que (1) seguir anadiendo interferencia, (2) mover la
%     variable de DECISION lo bastante para cambiar el veredicto, y (3) hacerlo por
%     FISICA y no por un artefacto del modelo de patron.
fprintf('\n(b) EJE DE SATURACION\n');
fprintf('    Sub-criterios: (1) aportacion marginal >= %.1f%% | (2) recorrido del KPI >= %.1f dB |\n', tolMarg, tolRange);
fprintf('                   (3) no dominado por el suelo de lobulos\n\n');
fprintf('  %7s %18s %14s %18s %14s\n','nRings','marg. r1 [%]','marg. r1[dB]','marg. ffr [%]','marg. ffr[dB]');
for p = 1:nP
    R = RES{p};
    fprintf('  %7d %18.3f %14.4f %18.3f %14.4f\n', R.nRings, ...
        R.scheme{1}.marg_pct, R.scheme{1}.marg_dB, ...
        R.scheme{2}.marg_pct, R.scheme{2}.marg_dB);
end

idx5 = find(nRingsList == 5, 1);
idx6 = find(nRingsList == 6, 1);
m5 = RES{idx5}.scheme{1};   m6 = RES{idx6}.scheme{1};
fprintf('\n  Aportacion marginal del anillo 5: %.3f%% del total (%.4f dB de subida de I_intra)\n', ...
    m5.marg_pct, m5.marg_dB);
fprintf('  Aportacion marginal del anillo 6: %.3f%% del total (%.4f dB de subida de I_intra)\n', ...
    m6.marg_pct, m6.marg_dB);

isMaterial = (m5.marg_pct >= tolMarg) || (m6.marg_pct >= tolMarg);

% --- Sub-criterio 2: RECORRIDO del KPI (un eje de saturacion tiene que mover la
%     variable de DECISION, no solo el presupuesto de interferencia) -----------
kpiRange = nan(1,nC);
for c = 1:nC
    v  = arrayfun(@(p) RES{p}.scheme{c}.SINR_edge_p5, 1:nP);
    bi = arrayfun(@(p) RES{p}.biased, 1:nP);
    vv = v(~bi);
    kpiRange(c) = max(vv) - min(vv);
    fprintf('\n  Recorrido de SINR_edge_p5 en %s (puntos NO sesgados): %.3f dB', ...
        cases{c}.name, kpiRange(c));
end
% Referencia: el eje de densidad orbital de la FASE A movio reuse1 unos 7 dB
% (-11.41 -> -4.41 dB entre T=66 y T=4000). Ese es el orden de magnitud que
% caracteriza a un eje con recorrido util.
fprintf('\n  (referencia FASE A: el eje de densidad orbital movio reuse1 ~7 dB)\n');
hasRange = max(kpiRange) >= tolRange;

% --- Sub-criterio 3: ¿es FISICA o es el SUELO DE LOBULOS? -------------------
%  beam_gain_dB aplica un suelo PLANO: Grel = max(Bessel, sidelobe_floor_dB). En
%  cuanto el Bessel cae por debajo del suelo, TODOS los haces de los anillos
%  siguientes quedan clavados en el mismo valor, y como el anillo r tiene 6r haces
%  su aportacion AGREGADA crece linealmente con r en vez de decaer. A partir de ese
%  anillo, lo que anade nRings es ARTEFACTO DEL MODELO DE PATRON, no interferencia
%  fisica. Se localiza el anillo de entrada del suelo con la geometria nominal
%  (usuario en el centro del cluster, offset r*s).
HPBW  = cfg.radio.beamwidth3dB_deg;
flr   = cfg.radio.sidelobe_floor_dB;
sdeg  = BLmax.spacing_deg;
floorOnset = NaN;
fprintf('\n  Patron por anillo (geometria nominal: usuario en el centro, offset = r*s):\n');
fprintf('  %6s %8s %12s %12s %14s %10s\n','anillo','haces','offset[deg]','Bessel[dB]','con suelo[dB]','aporta rel');
for r = 0:max(nRingsList)
    nb  = max(6*r, 1);
    gb  = beam_gain_dB(r*sdeg, HPBW, -Inf);
    gf  = beam_gain_dB(r*sdeg, HPBW, flr);
    atF = gb < flr;
    if atF && isnan(floorOnset), floorOnset = r; end
    if atF, mk = '<- SUELO'; else, mk = ''; end
    fprintf('  %6d %8d %12.4f %12.2f %14.2f %10.4f  %s\n', r, nb, r*sdeg, gb, gf, nb*10^(gf/10), mk);
end
isFloorDriven = ~isnan(floorOnset) && floorOnset <= max(nRingsList);
if isFloorDriven
    fprintf(['\n  -> El suelo de lobulos (%.0f dB) entra en el anillo %d: de ahi en adelante\n' ...
             '     cada anillo aporta 6r haces CLAVADOS en el suelo y su aportacion agregada\n' ...
             '     CRECE con r en vez de decaer. Con el Bessel real esos anillos no aportarian\n' ...
             '     nada (%.1f dB en el anillo %d). La cola es ARTEFACTO DEL MODELO DE PATRON.\n'], ...
        flr, floorOnset, beam_gain_dB(max(nRingsList)*sdeg, HPBW, -Inf), max(nRingsList));
end

% --- VEREDICTO (b): los tres sub-criterios ---------------------------------
isAxis = isMaterial && hasRange && ~isFloorDriven;
fprintf('\n  Sub-criterios:\n');
fprintf('    1) aportacion marginal >= %.1f%% en n=5/6 ....... %s\n', tolMarg, yesno(isMaterial));
fprintf('    2) recorrido del KPI  >= %.1f dB ................ %s (max %.3f dB)\n', ...
    tolRange, yesno(hasRange), max(kpiRange));
fprintf('    3) NO dominado por el suelo de lobulos ......... %s\n', yesno(~isFloorDriven));

fprintf('\n  CONCLUSION (b): ');
if isAxis
    fprintf(['nRings escala la interferencia dominante con recorrido util y sin\n' ...
             '                  artefacto -> RECOMENDACION: USAR nRings COMO BARRIDO SECUNDARIO DE E3.\n']);
else
    fprintf('nRings NO sirve como eje de saturacion.\n');
    if ~hasRange
        fprintf(['                  Aunque la aportacion marginal no sea nula, el KPI solo se mueve\n' ...
                 '                  %.2f dB en TODO el rango util (x%.1f haces co-canal): el eje no\n' ...
                 '                  cambia el veredicto de viabilidad.\n'], ...
            max(kpiRange), RES{end}.nCo_reuse1 / RES{3}.nCo_reuse1);
    end
    if isFloorDriven
        fprintf(['                  Y lo poco que se mueve viene del SUELO DE LOBULOS PLANO, no de\n' ...
                 '                  fisica: saturar por ahi seria saturar un artefacto del modelo.\n']);
    end
    fprintf('                  RECOMENDACION: USAR nRings SOLO COMO VALIDACION.\n');
    fprintf(['                  El eje que SI escala la interferencia intra es la DENSIDAD DE HACES\n' ...
             '                  por unidad de area (tamano de celda: cfg.radio.beamwidth3dB_deg /\n' ...
             '                  separacion inter-haz) o cfg.geom.minElev (compresion angular), no\n' ...
             '                  anadir anillos con la separacion FIJA.\n']);
end
isNegligible = ~isMaterial;   % se conserva por compatibilidad del .mat

%% ------------------------------------------------------------------------
%  6. TAREA 3 - COSTE
%  ------------------------------------------------------------------------
fprintf('\n\n=======================================================================\n');
fprintf('  TAREA 3: COSTE POR PUNTO vs nRings\n');
fprintf('=======================================================================\n\n');
fprintf('(el tiempo NO incluye la tabla P.618: se construye antes del barrido)\n');
fprintf('%7s %8s %10s %16s %18s %16s\n', ...
    'nRings','nBeams','tiempo[s]','pico bloque[GB]','pico ventana[GB]','delta RAM[GB]');
fprintf('%s\n', repmat('-',1,80));
for p = 1:nP
    R = RES{p};
    fprintf('%7d %8d %10.1f %16.3f %18.3f %16.2f\n', ...
        R.nRings, R.nBeams, R.elapsed_s, R.peak_model_GB, R.peak_full_GB, R.peak_sampled_GB);
end
t1 = RES{1}.elapsed_s;  tN = RES{end}.elapsed_s;
fprintf('\nEl coste sube poco con nRings: el termino dominante del modelo es M*N*Ntb (geometria,\n');
fprintf('independiente del cluster) y el de haces es M*nBeams*Ntb. Tiempo x%.2f de nRings=%d a %d.\n', ...
    tN/max(t1,eps), nRingsList(1), nRingsList(end));

%% ------------------------------------------------------------------------
%  7. Guardado
%  ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%  6b. FIGURAS
%  Esta seccion no tenia ninguna, y es la que peor se entiende sin ver el
%  perfil de interferencia por anillo: el hallazgo (que la cola la fija el SUELO
%  DE LOBULOS PLANO y no interferencia fisica) es dificil de defender solo con
%  numeros. Se exportan a 300 dpi con save_fig, como el resto del proyecto.
%  ------------------------------------------------------------------------
FIGDIR = 'figs_nrings';

nRv = arrayfun(@(p) RES{p}.nRings, 1:nP);
bi  = arrayfun(@(p) RES{p}.biased, 1:nP);

% (a) CONVERGENCIA del KPI frente al tamano del cluster, con los puntos SESGADOS
%     marcados aparte: incluirlos en la lectura de convergencia seria un error.
fA = figure('Name','nRings | convergencia del KPI','Color','w');
mk = {'o-','s-'};
for c = 1:nC
    v = arrayfun(@(p) RES{p}.scheme{c}.SINR_edge_p5, 1:nP);
    plot(nRv(~bi), v(~bi), mk{min(c,2)}, 'LineWidth',1.7, ...
         'DisplayName',cases{c}.name); hold on;
    if any(bi)
        plot(nRv(bi), v(bi), 'x', 'MarkerSize',12, 'LineWidth',2, ...
             'Color',[.7 .1 .1], 'HandleVisibility','off');
    end
end
if any(bi)
    plot(nan, nan, 'x', 'MarkerSize',12, 'LineWidth',2, 'Color',[.7 .1 .1], ...
         'DisplayName','punto SESGADO (rejilla fuera del cluster)');
end
xline(5, 'k--', 'LineWidth',1.4, 'DisplayName','valor de trabajo (nRings=5)');
grid on; xlabel('nRings (anillos completos del cluster)');
ylabel('SINR_{edge,p5} [dB]');
title({'Convergencia del KPI de borde frente al tamano del cluster', ...
       'la curva se aplana, pero el residuo lo fija el SUELO DE LOBULOS (ver fig. c)'});
legend('Location','best');
save_fig(fA, FIGDIR, 'nrings_a_convergencia_kpi_vs_tamano_cluster');

% (b) PERFIL DE INTERFERENCIA POR ANILLO en el cluster mas grande: donde vive
%     realmente la interferencia intra-satelite.
fB = figure('Name','nRings | perfil por anillo','Color','w');
pr = RES{end}.scheme{1}.pct_ring;                    % reuse1, cluster mayor
bar(0:numel(pr)-1, pr, 'FaceColor',[0.30 0.55 0.75]); hold on;
cum = cumsum(pr);
plot(0:numel(pr)-1, cum, 'k.-', 'LineWidth',1.6, 'MarkerSize',14);
yline(100,'k:');
grid on; xlabel('anillo r (0 = haz propio)');
ylabel('aportacion a la interferencia INTRA [%]');
title(sprintf(['Perfil de interferencia por anillo (reuso-1, nRings=%d): ' ...
               'los anillos 0-4 son el %.1f%%'], RES{end}.nRings, cum(min(5,end))));
legend({'aportacion del anillo','acumulada'}, 'Location','east');
save_fig(fB, FIGDIR, 'nrings_b_perfil_interferencia_por_anillo');

% (c) LA FIGURA DEL HALLAZGO: patron Bessel REAL frente al SUELO PLANO, en la
%     posicion angular de cada anillo. Donde las dos curvas se cruzan es donde el
%     suelo empieza a MANDAR, y a partir de ahi la cola de interferencia es un
%     ARTEFACTO DEL MODELO, no fisica.
fC = figure('Name','nRings | el suelo de lobulos manda a partir del anillo 5','Color','w');
sdeg  = RES{end}.spacing_deg;
HPBWd = cfg.radio.beamwidth3dB_deg;
flr   = cfg.radio.sidelobe_floor_dB;
rr    = 0:RES{end}.nRings;
gBes  = arrayfun(@(r) beam_gain_dB(r*sdeg, HPBWd, -Inf), rr);   % Bessel PURO
gUsed = arrayfun(@(r) beam_gain_dB(r*sdeg, HPBWd, flr),  rr);   % con suelo
plot(rr, gBes,  'o-', 'LineWidth',1.8, 'DisplayName','Bessel real (sin suelo)'); hold on;
plot(rr, gUsed, 's--','LineWidth',1.8, 'DisplayName','patron USADO (con suelo)');
yline(flr, 'k:', 'LineWidth',1.6, ...
      'Label',sprintf('suelo de lobulos = %g dB', flr));
iOn = find(gBes < flr, 1);
if ~isempty(iOn)
    xline(rr(iOn), 'r--', 'LineWidth',1.8, ...
        'Label',sprintf('el suelo MANDA desde el anillo %d', rr(iOn)), ...
        'LabelVerticalAlignment','bottom');
end
grid on; xlabel('anillo r'); ylabel('ganancia relativa del haz co-canal [dB]');
title({'Por que nRings NO es un eje de saturacion', ...
       'a partir del cruce, los 6r haces de cada anillo quedan CLAVADOS en el suelo'});
legend('Location','southwest');
save_fig(fC, FIGDIR, 'nrings_c_suelo_de_lobulos_vs_bessel_real');

fprintf('\n[figuras] 3 PNG a 300 dpi en %s%s\n', FIGDIR, filesep);

timing = struct('warmup_p618_s', tWarm, 'study_s', toc(tStudy));
fprintf('\n[tiempos] estudio %.1f s (%.1f min) (+ %.1f s de warmup P.618)\n', ...
    timing.study_s, timing.study_s/60, timing.warmup_p618_s);

save('convergence_nrings.mat', 'RES', 'nRingsList', 'cases', 'cfg', 'timeBlock', ...
     'nmin', 'nminGlobal', 'isNegligible', 'isAxis', 'isMaterial', 'hasRange', ...
     'isFloorDriven', 'floorOnset', 'kpiRange', 'tolConv', 'tolMarg', 'tolRange', ...
     'timing', '-v7.3');
fprintf('\nResultados en convergence_nrings.mat\n');

%% ========================================================================
%  FUNCIONES LOCALES
%  ========================================================================
function s = yesno(b)
%YESNO  Etiqueta SI/NO para las tablas de veredicto.
if b, s = 'SI'; else, s = 'NO'; end
end

%  ------------------------------------------------------------------------
%  Descomposicion por anillo de la interferencia intra
%  ------------------------------------------------------------------------
function [Pr, Ptot, nCov, relErr] = ring_decomposition(ctx, A, F, BL)
%RING_DECOMPOSITION  Reparte la interferencia INTRA por anillo hexagonal.
%   [Pr, Ptot, nCov, relErr] = ring_decomposition(ctx, A, F, BL)
%
%   NO es un modelo nuevo: es la MISMA suma co-canal de compute_sinr_ffr
%       I_intra_lin(m,k) = sum_{b in co-canal(m,k)} Glin(m,b,k) - Glin(m,b_propio,k)
%   partida segun BL.ring, respetando la regla de co-canal del motor:
%       A.subband == 0  -> co-canal = TODOS los haces
%       A.subband == c  -> co-canal = solo los haces de color c
%   El haz propio se descuenta del anillo al que pertenece (es la portadora).
%
%   Devuelve potencias ABSOLUTAS acumuladas (suma sobre las muestras con cobertura):
%     Pr    [1 x nRings+1] potencia interferente lineal aportada por cada anillo [W]
%     Ptot  potencia intra total segun el MOTOR (10.^(F.I_intra_dBW/10)), para control
%     nCov  numero de muestras con cobertura sumadas
%     relErr  max error relativo |sum(Pr) - Ptot| / Ptot  -> valida la particion

Glin  = ctx.Glin;                          % [M x nB x Nt]
cov   = ctx.cov;
ringv = BL.ring(:);                        % [nB x 1] anillo de cada haz
nRmax = max(ringv);
color = A.color(:);
sb    = A.subband;                         % [M x Nt]
Delta = A.Delta;
[Mloc, nB, Ntloc] = size(Glin);

% Factor de potencia comun a portadora e interferencia intra (PSD * ancho):
%   I_dBW = Cpsd + 10*log10(Bu/1e6) + 10*log10(I_lin_rel)
pw = 10.^((ctx.Cpsd_dBWMHz + 10*log10(A.B_user_Hz/1e6))/10);   % [M x Nt]

% Ganancia del haz PROPIO y su anillo
[mI, kI] = ndgrid(1:Mloc, 1:Ntloc);
ownLin = nan(Mloc, Ntloc);
ownRing= nan(Mloc, Ntloc);
i3 = sub2ind([Mloc nB Ntloc], mI(cov), ctx.beam(cov), kI(cov));
ownLin(cov)  = Glin(i3);
ownRing(cov) = ringv(ctx.beam(cov));

Pr   = zeros(1, nRmax+1);
Psum = zeros(Mloc, Ntloc);

for r = 0:nRmax
    selR = (ringv == r);

    % Suma de ganancia del anillo r restringida al conjunto co-canal del usuario
    sR = nan(Mloc, Ntloc);

    % (i) sub-banda interior (subband = 0): co-canal = TODOS los haces del anillo
    isInt = cov & (sb == 0);
    if any(isInt(:))
        sAll = reshape(sum(Glin(:, selR, :), 2, 'omitnan'), Mloc, Ntloc);
        sR(isInt) = sAll(isInt);
    end
    % (ii) sub-banda de color c: solo los haces del anillo r con ese color
    for c = 1:Delta
        selC = cov & (sb == c);
        if ~any(selC(:)), continue; end
        sCol = reshape(sum(Glin(:, selR & (color == c), :), 2, 'omitnan'), Mloc, Ntloc);
        sR(selC) = sCol(selC);
    end

    % Descontar el haz propio del anillo que lo contiene (es la portadora, no I)
    isOwn = cov & (ownRing == r);
    sR(isOwn) = sR(isOwn) - ownLin(isOwn);
    sR = max(sR, 0);

    Pj      = pw .* sR;                     % potencia absoluta del anillo r [W]
    Pj(~cov)= 0;  Pj(isnan(Pj)) = 0;
    Pr(r+1) = sum(Pj(:));
    Psum    = Psum + Pj;
end

% Control contra el MOTOR: la suma de anillos debe reproducir I_intra
Imot = 10.^(F.I_intra_dBW/10);
Imot(~cov) = 0;  Imot(isnan(Imot)) = 0;
Ptot = sum(Imot(:));
nCov = sum(cov(:));

den    = max(Ptot, realmin);
relErr = abs(sum(Pr) - Ptot) / den;
end
