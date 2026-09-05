%% RUN_CONVERGENCE_STUDYA  FASE A: convergencia numerica y memoria (previo a E3).
%  NO es el barrido de saturacion. Es el estudio que fija los PARAMETROS con los
%  que se lanzara E3: cuanta rejilla de usuarios y cuantos instantes hacen falta
%  para que el percentil 5 de SINR de BORDE (la metrica reina de H1) sea ESTABLE,
%  y cuanta RAM cuesta cada realizacion.
%
%  ESCENARIO FIJO: densidad T = 1584 (P = 72), Starlink-like, INTRA-constelacion,
%  banda Ku (config_default). Solo se barren los parametros NUMERICOS del
%  experimento, nunca la fisica.
%
%  TAREAS
%    1. Convergencia de SINR_edge_p5 frente a la resolucion de la rejilla (step_km).
%    2. Convergencia de SINR_edge_p5 frente al numero de instantes (Nt).
%    3. Estabilidad entre REALIZACIONES (dispersion del p5).
%    4. Saturacion del numero de interferentes co-canal VISIBLES frente a T.
%    5. Informe de memoria y veredicto de viabilidad de E3 en esta maquina.
%
%  TODO EN SECUENCIAL (sin parfor) a proposito: esta fase mide el consumo de UN
%  solo proceso, que es justo el dato que despues decide cuantos workers caben.
%
%  ---------------------------------------------------------------------------
%  MODELO DE MEMORIA (el numero clave de esta fase)
%  ---------------------------------------------------------------------------
%  El pico NO son los tensores de compute_geometry. El maximo se alcanza DENTRO de
%  compute_interference, en la llamada  atm_loss_dB(cfg, G.el)  sobre un array
%  [M x N x Nt]. En ese instante coexisten, por ELEMENTO M*N*Nt:
%       G.el, G.az, G.range        3 x 8 = 24 B
%       G.vis (logical)                    1 B
%       FSPL_all                           8 B
%       Latm_clear, Latm_rain      2 x 8 = 16 B   (dentro de atm_loss_dB)
%       valid (logical)                    1 B
%       temporales de interp1      ~3 x 8 = 24 B  (el_deg(valid), interp1, max)
%       -----------------------------------------
%       PICO ~ 74 B por elemento M*N*Nt
%  A eso se suma ffr_context, que es M*nBeams*Nt*(psi+Grel+Glin) ~ 40 B/elem pero
%  con nBeams (91) en lugar de N (1584): en este escenario aporta < 4% y no manda.
%  Con 553 usuarios, T=1584 y Nt=61 el pico modelado es ~3.9 GB; con Nt=241, ~15.6 GB.
%
%  ---------------------------------------------------------------------------
%  TROCEADO TEMPORAL (por que es EXACTO y no una aproximacion)
%  ---------------------------------------------------------------------------
%  Varias celdas del barrido pedido no caben de una pieza en 7.4 GB. En lugar de
%  reducir el experimento, este script evalua el pipeline por BLOQUES DE INSTANTES
%  y concatena las matrices [M x Nt] antes de llamar a compute_kpis. Es exacto
%  porque NADA del pipeline acopla instantes distintos:
%     - geometria, radioenlace e interferencia se calculan por instante;
%     - ffr_policy fija tau por INSTANTE (modo cuantil) y alpha por INSTANTE;
%     - ffr_allocate mide la carga n0 por INSTANTE;
%     - compute_kpis toma percentiles sobre la UNION de muestras (m,k) y la unica
%       agregacion temporal (R_agg) es una media sobre instantes.
%  El troceado se VALIDA al arrancar (seccion 1): mismo caso troceado y sin trocear
%  debe dar max|dif| = 0. Si no lo diera, el script para.
%  NO se toca el motor: el troceado vive aqui, en el script de estudio.
%
%  ---------------------------------------------------------------------------
%  SEGURIDAD DE MEMORIA
%  ---------------------------------------------------------------------------
%  Antes de cada caso se estima el pico y se elige el tamano de bloque para no
%  superar ST.mem_frac de la RAM fisica LIBRE. Si ni siquiera UN instante cabe, el
%  caso se SALTA con un mensaje explicito y el estudio continua: nunca se cuelga la
%  maquina ni se aborta el estudio entero.
%
%  Uso: ejecutar en la raiz del proyecto. Resultados en convergence_studyA.mat.
%
%  Aportacion propia (metodologia experimental; no altera la fisica).

clear; clc; close all;

%% 0. Parametros del estudio                                    % <-- DECISION
ST.T_fixed   = 1584;      % densidad fija de las tareas 1-3 (Starlink S1)
ST.P_fixed   = 72;
ST.radius_km = 40;        % radio de la rejilla local (fijo en toda la tarea 1)
ST.window_s  = 7200;      % ventana temporal comun (2 h)

ST.steps_km  = [8 6 4 3 2];          % tarea 1: resolucion de la rejilla
ST.Nt_task1  = 61;                   % tarea 1: instantes (dt = 120 s)

ST.Nt_list   = [31 61 121 241];      % tarea 2: instantes (misma ventana, otro dt)

ST.nReal     = 5;                    % tarea 3: realizaciones (fase orbital)

ST.T_list    = [66 300 600 1000 1584 2500 4000];   % tarea 4: densidades
ST.P_list    = [ 3  12  24   40   72  100  160];   % P ~ T/22..25 (T multiplo de P)
ST.Nt_task4  = 5;                    % tarea 4: instantes (media estable de visibles)

ST.tol_dB     = 0.2;      % criterio de convergencia del p5
ST.tol_std_dB = 0.3;      % criterio de fiabilidad (dispersion entre realizaciones)

ST.mem_frac        = 0.55;   % fraccion de la RAM LIBRE que se permite ocupar
ST.bytes_per_elem  = 74;     % modelo del pico por elemento M*N*Nt (ver cabecera)
ST.bytes_per_belem = 40;     % idem por elemento M*nBeams*Nt (ffr_context)
ST.elMin_valid     = 45;     % deg: submuestra de GEOMETRIA VALIDA (informativa)

%% 0b. Configuracion base (perfil de estudio Ku, identico a run_ffr_demo)
cfg0 = config_default();
cfg0.ground.mode      = 'localgrid';
cfg0.ground.radius_km = ST.radius_km;
cfg0.ground.step_km   = 3;
cfg0.time.t0          = 0;
cfg0.time.duration    = ST.window_s;
cfg0.ffr.tau_mode     = 'quantile';
cfg0.ffr.tau_q        = 50;
cfg0.ffr.sched        = 'share';
cfg0.beams.nRings     = 5;                 % 91 haces
cfg0.constellations(1).T = ST.T_fixed;
cfg0.constellations(1).P = ST.P_fixed;

% Los dos esquemas que se comparan en las tareas 1-3
SCH = { struct('scheme','reuse1','Delta',1,'alpha',NaN,'name','reuse1'), ...
        struct('scheme','ffr',   'Delta',3,'alpha',0.4,'name','ffr(D=3,a=0.4)') };

fprintf('=========================================================================\n');
fprintf('  FASE A: ESTUDIO DE CONVERGENCIA Y MEMORIA  (paso previo a E3)\n');
fprintf('=========================================================================\n');
fprintf('Escenario fijo: T=%d (P=%d) intra-constelacion Ku | radio %d km | ventana %.0f h\n', ...
    ST.T_fixed, ST.P_fixed, ST.radius_km, ST.window_s/3600);
fprintf('Ejecucion SECUENCIAL (sin parfor): se mide el consumo de UN proceso.\n');
fprintf('RAM fisica libre al arrancar: %.2f GB | presupuesto (%.0f%%): %.2f GB\n', ...
    ram_free_GB(), 100*ST.mem_frac, ST.mem_frac*ram_free_GB());

% Precalentar la tabla P.618 fuera de todo cronometro (coste fijo ~80 s por proceso)
fprintf('\n[warmup] construyendo la tabla ITU-R P.618...\n');
tw = tic;  atm_loss_dB(cfg0, 45);  tWarm = toc(tw);
fprintf('[warmup] listo (%.1f s)\n', tWarm);

% CRONOMETRO DE EXTREMO A EXTREMO del estudio, arrancado
% DESPUES del warmup: la tabla P.618 es un coste fijo de la configuracion
% atmosferica, no del estudio, y se reporta aparte.
tStudy = tic;

% Linea base de memoria del proceso: todo lo que se mida despues es lo que
% realmente cuesta el pipeline, sin contar MATLAB ni la tabla P.618.
pk_base = proc_peak_GB();
fprintf('[memoria] linea base del proceso tras el warmup: %.2f GB\n', pk_base);

%% 1. VALIDACION DEL TROCEADO TEMPORAL (debe ser EXACTO)
%  Caso pequeno que cabe entero: se evalua sin trocear y troceado en bloques de 2
%  instantes. Si el troceado no fuese exacto, todo el estudio quedaria invalidado.
fprintf('\n----- [0] Validacion del troceado temporal (exactitud) -----\n');
cfgV = cfg0;
cfgV.ground.step_km = 8;
cfgV.constellations(1).T = 300;  cfgV.constellations(1).P = 12;
cfgV.time.duration = 1800;  cfgV.time.dt = 300;      % Nt = 7

Efull = eval_case(cfgV, SCH, ST, Inf);      % sin trocear
Echnk = eval_case(cfgV, SCH, ST, 2);        % bloques de 2 instantes
[mdC, whC] = deep_maxdiff({Efull.K}, {Echnk.K}, 'K');
fprintf('bloques: %d (sin trocear) vs %d (troceado) | max|dif| = %.3g  (%s)\n', ...
    Efull.nChunks, Echnk.nChunks, mdC, whC);
if ~(mdC == 0)
    error('run_convergence_studyA:chunk', ...
        ['El troceado temporal NO es exacto (max|dif| = %.6g en %s). ' ...
         'PARAR: el resto del estudio no seria fiable.'], mdC, whC);
end
fprintf('OK: el troceado es EXACTO (max|dif| = 0). Se puede usar para acotar memoria.\n');

%% 2. TAREA 1 -- Convergencia frente al tamano de rejilla
fprintf('\n=========== TAREA 1: convergencia vs rejilla (step_km) ===========\n');
fprintf('T=%d, Nt=%d (dt=%.0f s), radio=%d km\n', ...
    ST.T_fixed, ST.Nt_task1, ST.window_s/(ST.Nt_task1-1), ST.radius_km);
fprintf('\n%-8s %7s %7s %10s %10s %10s %10s %9s %9s %7s\n', ...
    'step_km','M','M_borde','p5_reuse1','p5_ffr','p5_r1(V)','p5_ffr(V)', ...
    'tens[GB]','pico[GB]','bloques');
fprintf('%s\n', repmat('-',1,100));

T1 = struct('step_km',ST.steps_km, 'M',nan(1,numel(ST.steps_km)), ...
            'Medge',nan(1,numel(ST.steps_km)), 'p5',nan(2,numel(ST.steps_km)), ...
            'p5V',nan(2,numel(ST.steps_km)), 'tens_GB',nan(1,numel(ST.steps_km)), ...
            'peak_model_GB',nan(1,numel(ST.steps_km)), 'peak_meas_GB',nan(1,numel(ST.steps_km)), ...
            'peak_abs_GB',nan(1,numel(ST.steps_km)), 'nChunks',nan(1,numel(ST.steps_km)), ...
            'ok',false(1,numel(ST.steps_km)), 'secs',nan(1,numel(ST.steps_km)));

for i = 1:numel(ST.steps_km)
    cfg = cfg0;
    cfg.ground.step_km = ST.steps_km(i);
    cfg.time.dt        = ST.window_s / (ST.Nt_task1 - 1);

    E = eval_case(cfg, SCH, ST);
    T1.M(i)             = E.M;
    T1.tens_GB(i)       = E.tens_GB;
    T1.peak_model_GB(i) = E.peak_model_GB;
    T1.peak_meas_GB(i)  = E.peak_meas_GB;
    T1.peak_abs_GB(i)   = E.peak_abs_GB;
    T1.nChunks(i)       = E.nChunks;
    T1.secs(i)          = E.secs;
    T1.ok(i)            = E.ok;

    if ~E.ok
        fprintf('%-8g %7d %7s %10s %10s %10s %10s %9.2f %9.2f %7s  << SALTADO: %s\n', ...
            ST.steps_km(i), E.M, '-','-','-','-','-', E.tens_GB, E.peak_model_GB, '-', E.why);
        continue;
    end

    T1.Medge(i)  = E.K{1}.viab.n_edge;
    T1.p5(:,i)   = [E.K{1}.viab.SINR_edge_p5;  E.K{2}.viab.SINR_edge_p5];
    T1.p5V(:,i)  = [E.KV{1}.viab.SINR_edge_p5; E.KV{2}.viab.SINR_edge_p5];
    fprintf('%-8g %7d %7d %10.3f %10.3f %10.3f %10.3f %9.2f %9.2f %7d\n', ...
        ST.steps_km(i), E.M, T1.Medge(i), T1.p5(1,i), T1.p5(2,i), ...
        T1.p5V(1,i), T1.p5V(2,i), E.tens_GB, E.peak_model_GB, E.nChunks);
end
fprintf('%s\n', repmat('-',1,100));

[step_conv, T1.dp5] = converged_step(ST.steps_km, T1.p5, ST.tol_dB, T1.ok);
fprintf('Diferencias |p5(step_i) - p5(step_{i+1})| (max de los dos esquemas) [dB]:');
for i = 1:numel(T1.dp5), fprintf(' %.3f', T1.dp5(i)); end
fprintf('\n');
if isnan(step_conv)
    fprintf('>> NO se alcanza el criterio de %.1f dB en el rango barrido.\n', ST.tol_dB);
else
    fprintf('>> step_km CONVERGIDO = %g km (al refinar mas, el p5 cambia < %.1f dB)\n', ...
        step_conv, ST.tol_dB);
end

%% 3. TAREA 2 -- Convergencia frente al numero de instantes
step_use = step_conv;  if isnan(step_use), step_use = min(ST.steps_km); end
fprintf('\n=========== TAREA 2: convergencia vs numero de instantes (Nt) ===========\n');
fprintf('T=%d, step_km=%g (convergido), ventana %.0f h fija -> cambia dt\n', ...
    ST.T_fixed, step_use, ST.window_s/3600);
fprintf('\n%-6s %8s %8s %10s %10s %10s %10s %9s %9s %7s\n', ...
    'Nt','dt[s]','M_borde','p5_reuse1','p5_ffr','p5_r1(V)','p5_ffr(V)','tens[GB]','pico[GB]','bloques');
fprintf('%s\n', repmat('-',1,96));

T2 = struct('Nt',ST.Nt_list, 'p5',nan(2,numel(ST.Nt_list)), 'p5V',nan(2,numel(ST.Nt_list)), ...
            'tens_GB',nan(1,numel(ST.Nt_list)), 'peak_model_GB',nan(1,numel(ST.Nt_list)), ...
            'nChunks',nan(1,numel(ST.Nt_list)), 'ok',false(1,numel(ST.Nt_list)), ...
            'secs',nan(1,numel(ST.Nt_list)), 'step_km',step_use);

for i = 1:numel(ST.Nt_list)
    cfg = cfg0;
    cfg.ground.step_km = step_use;
    cfg.time.dt        = ST.window_s / (ST.Nt_list(i) - 1);

    E = eval_case(cfg, SCH, ST);
    T2.tens_GB(i)       = E.tens_GB;
    T2.peak_model_GB(i) = E.peak_model_GB;
    T2.nChunks(i)       = E.nChunks;
    T2.secs(i)          = E.secs;
    T2.ok(i)            = E.ok;
    if ~E.ok
        fprintf('%-6d %8.0f %8s %10s %10s %10s %10s %9.2f %9.2f %7s  << SALTADO: %s\n', ...
            ST.Nt_list(i), cfg.time.dt, '-','-','-','-','-', E.tens_GB, E.peak_model_GB, '-', E.why);
        continue;
    end
    T2.p5(:,i)  = [E.K{1}.viab.SINR_edge_p5;  E.K{2}.viab.SINR_edge_p5];
    T2.p5V(:,i) = [E.KV{1}.viab.SINR_edge_p5; E.KV{2}.viab.SINR_edge_p5];
    fprintf('%-6d %8.0f %8d %10.3f %10.3f %10.3f %10.3f %9.2f %9.2f %7d\n', ...
        ST.Nt_list(i), cfg.time.dt, E.K{1}.viab.n_edge, T2.p5(1,i), T2.p5(2,i), ...
        T2.p5V(1,i), T2.p5V(2,i), E.tens_GB, E.peak_model_GB, E.nChunks);
end
fprintf('%s\n', repmat('-',1,96));

[Nt_conv, T2.dp5] = converged_step(ST.Nt_list, T2.p5, ST.tol_dB, T2.ok);
fprintf('Diferencias |p5(Nt_i) - p5(Nt_{i+1})| (max de los dos esquemas) [dB]:');
for i = 1:numel(T2.dp5), fprintf(' %.3f', T2.dp5(i)); end
fprintf('\n');
if isnan(Nt_conv)
    fprintf('>> NO se alcanza el criterio de %.1f dB en el rango barrido.\n', ST.tol_dB);
else
    fprintf('>> Nt CONVERGIDO = %d instantes (al doblarlo, el p5 cambia < %.1f dB)\n', ...
        Nt_conv, ST.tol_dB);
end

%% 4. TAREA 3 -- Estabilidad entre realizaciones
Nt_use = Nt_conv;  if isnan(Nt_use), Nt_use = ST.Nt_task1; end
fprintf('\n=========== TAREA 3: estabilidad entre realizaciones ===========\n');
fprintf('step_km=%g, Nt=%d, T=%d. La REALIZACION se cambia moviendo la FASE ORBITAL\n', ...
    step_use, Nt_use, ST.T_fixed);
fprintf(['(M0 y Om0 de la constelacion), NO el punto de tierra: mover cfg.ground.point\n' ...
         'cambiaria la clave de la cache P.618 y reconstruiria la tabla en cada realizacion.\n']);

% Realizacion 1 = nominal; el resto, fases sorteadas con semilla FIJA (reproducible
% y no elegida a mano).
rs = RandStream('Threefry','Seed',1);
M0v  = [0, 360*rand(rs,1,ST.nReal-1)];
Om0v = [0, 360*rand(rs,1,ST.nReal-1)];

T3 = struct('M0',M0v, 'Om0',Om0v, 'p5',nan(2,ST.nReal), 'p5V',nan(2,ST.nReal), ...
            'ok',false(1,ST.nReal), 'step_km',step_use, 'Nt',Nt_use);

fprintf('\n%-6s %9s %9s %12s %12s %12s %12s\n', ...
    'real.','M0[deg]','Om0[deg]','p5_reuse1','p5_ffr','p5_r1(V)','p5_ffr(V)');
fprintf('%s\n', repmat('-',1,76));
for i = 1:ST.nReal
    cfg = cfg0;
    cfg.ground.step_km = step_use;
    cfg.time.dt        = ST.window_s / (Nt_use - 1);
    cfg.constellations(1).M0  = M0v(i);
    cfg.constellations(1).Om0 = Om0v(i);

    E = eval_case(cfg, SCH, ST);
    T3.ok(i) = E.ok;
    if ~E.ok
        fprintf('%-6d %9.2f %9.2f %12s %12s %12s %12s  << SALTADO: %s\n', ...
            i, M0v(i), Om0v(i), '-','-','-','-', E.why);
        continue;
    end
    T3.p5(:,i)  = [E.K{1}.viab.SINR_edge_p5;  E.K{2}.viab.SINR_edge_p5];
    T3.p5V(:,i) = [E.KV{1}.viab.SINR_edge_p5; E.KV{2}.viab.SINR_edge_p5];
    fprintf('%-6d %9.2f %9.2f %12.3f %12.3f %12.3f %12.3f\n', ...
        i, M0v(i), Om0v(i), T3.p5(1,i), T3.p5(2,i), T3.p5V(1,i), T3.p5V(2,i));
end
fprintf('%s\n', repmat('-',1,76));

T3.mean = mean(T3.p5(:,T3.ok), 2);
T3.std  = std( T3.p5(:,T3.ok), 0, 2);
T3.meanV= mean(T3.p5V(:,T3.ok), 2);
T3.stdV = std( T3.p5V(:,T3.ok), 0, 2);
fprintf('%-18s  media = %7.3f dB   desv.tipica = %6.3f dB\n', 'reuse1',        T3.mean(1), T3.std(1));
fprintf('%-18s  media = %7.3f dB   desv.tipica = %6.3f dB\n', 'ffr(D=3)',      T3.mean(2), T3.std(2));
fprintf('%-18s  media = %7.3f dB   desv.tipica = %6.3f dB\n', 'reuse1 (V)',    T3.meanV(1),T3.stdV(1));
fprintf('%-18s  media = %7.3f dB   desv.tipica = %6.3f dB\n', 'ffr(D=3) (V)',  T3.meanV(2),T3.stdV(2));
T3.reliable = all(T3.std < ST.tol_std_dB);
fprintf('>> Criterio de fiabilidad (desv. tipica < %.1f dB): %s\n', ...
    ST.tol_std_dB, ternary_str(T3.reliable, 'CUMPLE', 'NO CUMPLE -> subir la rejilla'));

%% 5. TAREA 4 -- Saturacion de interferentes co-canal visibles vs densidad
Nt4 = Nt_use;      % se usa el Nt CONVERGIDO, no un solo instante
fprintf('\n=========== TAREA 4: interferentes co-canal VISIBLES vs densidad ===========\n');
fprintf(['step_km=%g, Nt=%d instantes (el CONVERGIDO de la tarea 2, no 1 solo:\n' ...
         'con 5 instantes el ruido de muestreo del p5 era ~0.4 dB, mayor que el efecto\n' ...
         'que se quiere medir; DECLARADO). P se escala con T (~22-25 sat/plano).\n'], ...
    step_use, Nt4);
fprintf('\n%-7s %6s %7s %11s %11s %10s %10s %11s %10s %8s\n', ...
    'T','P','cobert.','nVis_medio','n_cocanal','p5_reuse1','p5_ffr', ...
    'I_inter/N','penal_r1','pico[GB]');
fprintf('%-7s %6s %7s %11s %11s %10s %10s %11s %10s %8s\n', ...
    '','','','','','[dB]','[dB]','[dB]','[dB]','');
fprintf('%s\n', repmat('-',1,100));

nT = numel(ST.T_list);
T4 = struct('T',ST.T_list, 'P',ST.P_list, 'nVis',nan(1,nT), 'nCo',nan(1,nT), ...
            'p5',nan(2,nT), 'cov',nan(1,nT), 'peak_model_GB',nan(1,nT), ...
            'INR',nan(1,nT), 'pen',nan(1,nT), ...
            'nChunks',nan(1,nT), 'ok',false(1,nT), 'step_km',step_use, 'Nt',Nt4);

for i = 1:nT
    cfg = cfg0;
    cfg.ground.step_km = step_use;
    cfg.constellations(1).T = ST.T_list(i);
    cfg.constellations(1).P = ST.P_list(i);
    cfg.time.dt = ST.window_s / (Nt4 - 1);

    E = eval_case(cfg, SCH, ST);
    T4.peak_model_GB(i) = E.peak_model_GB;
    T4.nChunks(i)       = E.nChunks;
    T4.ok(i)            = E.ok;
    if ~E.ok
        fprintf('%-7d %6d %7s %11s %11s %10s %10s %11s %10s %8.2f  << SALTADO: %s\n', ...
            ST.T_list(i), ST.P_list(i), '-','-','-','-','-','-','-', E.peak_model_GB, E.why);
        continue;
    end
    T4.nVis(i) = E.nVis_mean;
    T4.nCo(i)  = E.nVis_mean - 1;              % co-canal = visibles menos el servidor
    T4.cov(i)  = E.covFrac;
    T4.p5(:,i) = [E.K{1}.viab.SINR_edge_p5; E.K{2}.viab.SINR_edge_p5];
    T4.INR(i)  = E.INR_dB(1);                  % I_inter / N en reuso-1
    T4.pen(i)  = E.K{1}.penalty_mean;          % penalizacion total C/N - SINR
    fprintf('%-7d %6d %7.2f %11.2f %11.2f %10.3f %10.3f %11.1f %10.2f %8.2f\n', ...
        ST.T_list(i), ST.P_list(i), T4.cov(i), T4.nVis(i), T4.nCo(i), ...
        T4.p5(1,i), T4.p5(2,i), T4.INR(i), T4.pen(i), E.peak_model_GB);
end
fprintf('%s\n', repmat('-',1,100));

% "Aplanamiento": el numero de visibles satura cuando el incremento RELATIVO al
% doblar (aprox.) la densidad cae por debajo del 10%.
T4.Tsat = NaN;
ok = T4.ok & ~isnan(T4.nCo);
Ts  = ST.T_list(ok);  Ns = T4.nCo(ok);
if numel(Ns) >= 2
    rel = diff(Ns) ./ max(Ns(1:end-1), eps);       % incremento relativo por tramo
    dens= diff(Ts) ./ Ts(1:end-1);                 % incremento relativo de densidad
    T4.elasticity = rel ./ max(dens, eps);         % d(nCo)/nCo por d(T)/T
    fprintf('Elasticidad d(n_cocanal)/n  por  d(T)/T  (1.0 = crece proporcional a T):');
    for i = 1:numel(T4.elasticity), fprintf(' %.2f', T4.elasticity(i)); end
    fprintf('\n');
    isat = find(T4.elasticity < 0.5, 1);           % crece menos de la mitad que T
    if ~isempty(isat), T4.Tsat = Ts(isat); end
end
if isnan(T4.Tsat)
    fprintf(['>> El RECUENTO de interferentes visibles NO se aplana: crece LINEALMENTE con T\n' ...
             '   (fraccion visible ~%.3f%% de la constelacion, constante). Es una identidad\n' ...
             '   geometrica -- el casquete de visibilidad sobre minElev=%g deg es un %% fijo\n' ...
             '   de la esfera -- , no un artefacto: por ese lado NO hay techo natural.\n'], ...
        100*mean(T4.nVis(ok)./ST.T_list(ok)), cfg0.geom.minElev);
else
    fprintf('>> Los interferentes visibles se APLANAN a partir de T ~ %d.\n', T4.Tsat);
end

% LA PREGUNTA OPERATIVA: no es cuantos interferentes se ven, sino a partir de que T
% deja de moverse el KPI. Es lo que fija el techo UTIL del barrido de E3.
[T4.Tsat_kpi, T4.dp5] = converged_step(ST.T_list, T4.p5, ST.tol_dB, T4.ok);
fprintf('\nDiferencias |p5(T_i) - p5(T_{i+1})| (max de los dos esquemas) [dB]:');
for i = 1:numel(T4.dp5), fprintf(' %.3f', T4.dp5(i)); end
fprintf('\n');
if isnan(T4.Tsat_kpi)
    fprintf('>> El KPI (SINR_edge_p5) NO se estabiliza en el rango barrido (tol %.1f dB).\n', ST.tol_dB);
else
    fprintf(['>> El KPI (SINR_edge_p5) se ESTABILIZA a partir de T ~ %d: mas alla, anadir\n' ...
             '   satelites NO cambia la SINR de borde.\n'], T4.Tsat_kpi);
end
fprintf(['>> I_inter/N maximo en todo el barrido: %.1f dB (T=%d). Si es << 0 dB, la\n' ...
         '   interferencia INTER-satelite sigue siendo DESPRECIABLE frente al ruido incluso\n' ...
         '   con %d satelites, y la dominante es la INTRA-satelite (fijada por el layout de\n' ...
         '   %d haces, que NO depende de T).\n'], ...
    max(T4.INR(ok)), ST.T_list(find(T4.INR == max(T4.INR(ok)), 1)), max(ST.T_list(ok)), ...
    build_beam_layout(cfg0).nBeams);

%% 6. TAREA 5 -- Informe de memoria y viabilidad de E3
fprintf('\n=========== TAREA 5: informe de memoria y viabilidad de E3 ===========\n');
freeGB  = ram_free_GB();
totalGB = ram_total_GB();

% Configuracion convergida (la que usara E3)
cfgE3 = cfg0;
cfgE3.ground.step_km = step_use;
cfgE3.time.dt        = ST.window_s / (Nt_use - 1);
uE3 = build_user_grid(cfgE3);
BLE3 = build_beam_layout(cfgE3);
mem1 = peak_model_GB(uE3.M, ST.T_fixed, Nt_use, BLE3.nBeams, ST);

fprintf('RAM total %.1f GB | libre ahora %.2f GB\n', totalGB, freeGB);
fprintf('\nConfiguracion CONVERGIDA para E3: step_km=%g (%d usuarios), Nt=%d, %d haces\n', ...
    step_use, uE3.M, Nt_use, BLE3.nBeams);
fprintf('%-46s %10s %10s\n','magnitud','[GB]','');
fprintf('%s\n', repmat('-',1,68));
fprintf('%-46s %10.2f\n', 'tensores G.el/az/range/vis (M*N*Nt*25 B)', ...
    uE3.M*ST.T_fixed*Nt_use*25/2^30);
fprintf('%-46s %10.2f\n', 'PICO modelado por realizacion (74 B/elem)', mem1);
fprintf('%-46s %10.2f\n', 'linea base del proceso (MATLAB + P.618)', pk_base);
fprintf('%-46s %10.2f\n', 'pico REAL del proceso en todo el estudio', proc_peak_GB());
fprintf('%s\n', repmat('-',1,68));

% VALIDACION EMPIRICA del modelo de 74 B/elemento: se usa el caso MAS GRANDE que
% se haya ejecutado SIN TROCEAR (nChunks == 1), que es el unico en el que el pico
% del proceso corresponde a una pasada entera. Si el modelo se quedara corto, el
% dimensionamiento de E3 seria optimista, asi que conviene contrastarlo.
iu = find(T1.ok & T1.nChunks == 1);
if ~isempty(iu)
    [~, jj] = max(T1.M(iu));  iu = iu(jj);
    elemU = T1.M(iu) * ST.T_fixed * ST.Nt_task1;
    bmeas = (T1.peak_abs_GB(iu) - pk_base) * 2^30 / elemU;
    fprintf(['[validacion del modelo] caso mas grande SIN trocear: step_km=%g ' ...
             '(%d usuarios, %d bloques)\n'], ST.steps_km(iu), T1.M(iu), T1.nChunks(iu));
    fprintf(['                        modelo %d B/elem  vs  medido %.0f B/elem ' ...
             '(pico %.2f GB sobre la linea base)\n'], ...
        ST.bytes_per_elem, bmeas, T1.peak_abs_GB(iu) - pk_base);
    if bmeas > ST.bytes_per_elem
        fprintf(['                        AVISO: el modelo se queda CORTO; usar %.0f B/elem\n' ...
                 '                        al dimensionar E3 (margen adicional x%.2f).\n'], ...
            bmeas, bmeas/ST.bytes_per_elem);
    else
        fprintf('                        el modelo es CONSERVADOR (cota superior): OK.\n');
    end
    memV.bytes_meas = bmeas;
else
    fprintf('[validacion del modelo] ningun caso se ejecuto sin trocear: sin validacion empirica.\n');
    memV.bytes_meas = NaN;
end

% --- A2: DIMENSIONADO DEL BARRIDO con el estimador COMPARTIDO ---------------
% estimate_sweep_memory es el que usan de verdad run_sweep_points / E3, y aplica
% ademas el troceado (cfg.compute.timeBlock). Se le pasa el punto MAS CARO del eje
% T para que el veredicto sea el del peor caso, y se contrasta con el pico que
% acaba de MEDIRSE en este mismo proceso: si el estimador y lo medido divergieran,
% el dimensionamiento de cualquier barrido futuro seria falso.
fprintf('\n--- Dimensionado del barrido (estimate_sweep_memory, fuente compartida) ---\n');
ptsE3 = arrayfun(@(t,p) struct('T',t,'P',p), ST.T_list, ST.P_list, 'UniformOutput', false);

cfgSweep = cfgE3;                       % SIN trocear (referencia)
cfgSweep.compute.timeBlock = [];
MEMfull = estimate_sweep_memory(cfgSweep, ptsE3, 4, struct('verbose',false));

fprintf('%-42s %10s %10s\n', 'timeBlock', 'pico/wk[GB]', 'caben');
fprintf('%s\n', repmat('-',1,64));
fprintf('%-42s %10.2f %10d\n', 'sin trocear (ventana entera)', ...
    MEMfull.peak_worker_GB, MEMfull.nFit);
tbList = [40 20 10 5];
MEMtb = cell(1,numel(tbList));
for i = 1:numel(tbList)
    c = cfgE3;  c.compute.timeBlock = tbList(i);
    MEMtb{i} = estimate_sweep_memory(c, ptsE3, 4, struct('verbose',false));
    fprintf('%-42s %10.2f %10d\n', sprintf('timeBlock = %d (%d bloques)', ...
        MEMtb{i}.nBlk, MEMtb{i}.nBlocks), MEMtb{i}.peak_worker_GB, MEMtb{i}.nFit);
end
fprintf('%s\n', repmat('-',1,64));

% timeBlock MINIMO que permite 4 workers, con la misma regla del estimador
i4 = find(cellfun(@(m) m.nFit >= 4, MEMtb), 1);
if isempty(i4)
    tb_reco = MEMfull.timeBlock_reco;
    fprintf('Para 4 workers hace falta timeBlock = %s (recomendacion del estimador).\n', ...
        ternary_str(isnumeric(tb_reco) && isscalar(tb_reco) && isfinite(tb_reco), ...
                    num2str(tb_reco), 'no alcanzable'));
else
    tb_reco = tbList(i4);
    fprintf('timeBlock RECOMENDADO para 4 workers: %d (pico/worker %.2f GB).\n', ...
        tb_reco, MEMtb{i4}.peak_worker_GB);
end
MEMA2 = struct('full',MEMfull, 'tbList',tbList, 'tb',{MEMtb}, 'tb_reco',tb_reco);

% CONTRASTE estimador vs MEDIDO, sobre el caso mas grande ejecutado sin trocear
if ~isnan(memV.bytes_meas)
    est_case = mem_peak_model_GB(T1.M(iu), ST.T_fixed, ST.Nt_task1, BLE3.nBeams);
    meas_case = T1.peak_abs_GB(iu) - pk_base;
    fprintf(['CONTRASTE estimador vs medido (caso sin trocear, %d usuarios): ' ...
             'modelo %.2f GB vs medido %.2f GB -> %+.1f%%\n'], ...
        T1.M(iu), est_case, meas_case, 100*(est_case-meas_case)/max(meas_case,eps));
    fprintf('  (el modelo debe quedar POR ENCIMA: es cota superior para dimensionar.)\n');
    MEMA2.est_case_GB = est_case;  MEMA2.meas_case_GB = meas_case;
end

nFit = floor(0.70*freeGB / max(mem1,eps));
fprintf('\nWorkers de parfor que caben (70%% de la RAM libre / pico por realizacion): %d\n', nFit);
if nFit >= 2
    verdict = sprintf('E3 VIABLE en esta maquina con %d worker(s)', min(nFit,4));
elseif nFit == 1
    verdict = 'E3 viable SOLO en secuencial (1 proceso); con mas RAM se paralelizaria';
else
    verdict = 'E3 NO cabe de una pieza: hace falta trocear en tiempo o mas RAM';
end
fprintf('VEREDICTO: %s\n', verdict);

% Que haria falta para 4 workers
need4 = 4*mem1/0.70;
fprintf('Para 4 workers en paralelo harian falta ~%.1f GB libres (%.1f GB de RAM total).\n', ...
    4*mem1, need4);

%% 7. RESUMEN FINAL
fprintf('\n=========================================================================\n');
fprintf('  RESUMEN FASE A  (parametros con los que se lanzara E3)\n');
fprintf('=========================================================================\n');
fprintf('  step_km recomendado ........ %s\n', ...
    ternary_str(isnan(step_conv), 'NO converge en el rango', sprintf('%g km  (%d usuarios)', step_use, uE3.M)));
fprintf('  Nt recomendado ............. %s\n', ...
    ternary_str(isnan(Nt_conv), 'NO converge en el rango', sprintf('%d instantes (dt = %.0f s)', Nt_use, ST.window_s/(Nt_use-1))));
fprintf('  desv. tipica del p5 ........ %.3f dB (reuse1) / %.3f dB (ffr)  -> %s\n', ...
    T3.std(1), T3.std(2), ternax(T3.reliable));
fprintf('  interferentes visibles ...... %s\n', ...
    ternary_str(isnan(T4.Tsat), 'NO saturan: crecen linealmente con T', sprintf('saturan en T ~ %d', T4.Tsat)));
fprintf('  saturacion del KPI .......... %s\n', ...
    ternary_str(isnan(T4.Tsat_kpi), 'el p5 no se estabiliza en el rango', ...
                sprintf('T ~ %d (mas alla, el p5 no se mueve)', T4.Tsat_kpi)));
fprintf('  I_inter/N maximo ............ %.1f dB (inter-satelite %s)\n', ...
    max(T4.INR(T4.ok)), ternary_str(max(T4.INR(T4.ok)) < -10, 'DESPRECIABLE', 'RELEVANTE'));
fprintf('  memoria por realizacion ..... %.2f GB (pico modelado)\n', mem1);
fprintf('  veredicto de memoria ........ %s\n', verdict);
fprintf('=========================================================================\n');

%% 8. Figuras
FIGDIR = 'figs_faseA';       % PNG a 300 dpi via save_fig (fuente unica)

fA = figure('Name','FASE A | convergencia vs rejilla','Color','w');
plot(ST.steps_km(T1.ok), T1.p5(1,T1.ok), 'o-','LineWidth',1.6,'DisplayName','reuse1'); hold on;
plot(ST.steps_km(T1.ok), T1.p5(2,T1.ok), 's-','LineWidth',1.6,'DisplayName','ffr (\Delta=3)');
set(gca,'XDir','reverse'); grid on;
xlabel('step\_km (rejilla mas FINA hacia la derecha)'); ylabel('SINR_{edge,p5} [dB]');
title('Convergencia del percentil 5 de borde frente a la resolucion de la rejilla');
legend('Location','best');
save_fig(fA, FIGDIR, 'faseA_a_convergencia_vs_paso_rejilla');

fB = figure('Name','FASE A | convergencia vs Nt','Color','w');
plot(ST.Nt_list(T2.ok), T2.p5(1,T2.ok), 'o-','LineWidth',1.6,'DisplayName','reuse1'); hold on;
plot(ST.Nt_list(T2.ok), T2.p5(2,T2.ok), 's-','LineWidth',1.6,'DisplayName','ffr (\Delta=3)');
grid on; xlabel('N_t (instantes)'); ylabel('SINR_{edge,p5} [dB]');
title('Convergencia del percentil 5 de borde frente al numero de instantes');
legend('Location','best');
save_fig(fB, FIGDIR, 'faseA_b_convergencia_vs_numero_instantes');

fC = figure('Name','FASE A | saturacion de interferentes','Color','w');
yyaxis left
plot(ST.T_list(T4.ok), T4.nCo(T4.ok), 'o-','LineWidth',1.6); ylabel('interferentes co-canal visibles');
yyaxis right
plot(ST.T_list(T4.ok), T4.p5(1,T4.ok), 's--','LineWidth',1.4); ylabel('SINR_{edge,p5} reuse1 [dB]');
grid on; xlabel('T (satelites)');
title('Saturacion del numero de interferentes co-canal visibles con la densidad');
save_fig(fC, FIGDIR, 'faseA_c_interferentes_y_kpi_vs_densidad_orbital');

fprintf('\n[figuras] 3 PNG a 300 dpi en %s%s\n', FIGDIR, filesep);

%% 9. Guardar
timing = struct('warmup_p618_s', tWarm, 'study_s', toc(tStudy));
fprintf('[tiempos] estudio %.1f s (%.1f min) (+ %.1f s de warmup P.618)\n', ...
    timing.study_s, timing.study_s/60, timing.warmup_p618_s);

save('convergence_studyA.mat','ST','cfg0','SCH','T1','T2','T3','T4', ...
     'step_conv','step_use','Nt_conv','Nt_use','mem1','nFit','verdict', ...
     'freeGB','totalGB','pk_base','memV','MEMA2','timing');
fprintf('\nResultados guardados en convergence_studyA.mat\n');

% =========================================================================
% ============================ FUNCIONES ==================================
% =========================================================================

function E = eval_case(cfg, SCH, ST, forceChunk)
%EVAL_CASE  Evalua un caso completo troceando el tiempo si hace falta.
%   Devuelve E.K / E.KV (KPIs por esquema), diagnostico y consumo de memoria.
%   E.ok = false si el caso no cabe ni con bloques de 1 instante (nunca cuelga).
if nargin < 4, forceChunk = []; end

tvec  = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
Nt    = numel(tvec);
sats  = build_constellation(cfg);
users = build_user_grid(cfg);
BL    = build_beam_layout(cfg);
M     = users.M;  N = sats.N;  nB = BL.nBeams;

E.M = M;  E.N = N;  E.Nt = Nt;  E.nBeams = nB;
E.tens_GB       = M*N*Nt*25 / 2^30;                       % tensores de geometria
E.peak_model_GB = peak_model_GB(M, N, Nt, nB, ST);        % pico de una pasada entera
E.ok = true;  E.why = '';

% --- Tamano de bloque temporal segun el presupuesto de RAM ---
budget_B = ST.mem_frac * ram_free_GB() * 2^30;
perInst  = M*(N*ST.bytes_per_elem + nB*ST.bytes_per_belem);   % bytes por instante
ntc = floor(budget_B / max(perInst, eps));
if ~isempty(forceChunk), ntc = min(ntc, forceChunk); end
if ntc < 1
    E.ok  = false;
    E.why = sprintf('ni 1 instante cabe (%.2f GB > %.2f GB disponibles)', ...
        perInst/2^30, budget_B/2^30);
    E.nChunks = NaN;  E.secs = 0;  E.peak_meas_GB = NaN;  E.peak_abs_GB = NaN;
    E.K = {};  E.KV = {};  E.nVis_mean = NaN;  E.covFrac = NaN;
    E.INR_dB = nan(1, numel(SCH));
    return;
end
ntc = min(ntc, Nt);
E.nChunks = ceil(Nt/ntc);

nS = numel(SCH);
% OJO CON EL TIPO de los acumuladores: en MATLAB  [ [] , logical ]  devuelve
% DOUBLE (el [] es double y double gana en la concatenacion). Si isCenter o cov
% degeneran a double, la indexacion logica posterior falla con "Array indices must
% be positive integers" en cuanto haya alguna muestra sin cobertura. Se inicializan
% por tanto como empties LOGICOS.
acc = repmat(struct('SINR',[],'CN',[],'pen',[],'Bu',[],'isC',logical([]),'alpha',[], ...
                    'Iint',[],'Nn',[]), 1, nS);
elevAll = [];  covAll = logical([]);  nVisAll = [];

pk0 = proc_peak_GB();
t0  = tic;

for c0 = 1:ntc:Nt
    idx = c0 : min(c0+ntc-1, Nt);
    tv  = tvec(idx);

    R_eci  = propagate(cfg, sats, tv);
    R_ecef = eci2ecef(cfg, R_eci, tv);
    % Mascara por constelacion (sats.minElev). Con una sola constelacion y
    % cfg.constellations(1).minElev = [] es IDENTICO a la mascara global, pero se
    % pasa explicitamente para no depender del camino por defecto.
    G      = compute_geometry(cfg, users, R_ecef, sats.minElev);
    S      = associate_serving(cfg, G);
    L      = compute_link_budget(cfg, S);
    INTF   = compute_interference(cfg, sats, users, G, S, L);
    ctx    = ffr_context(cfg, sats, users, BL, R_ecef, G, S, L, INTF, false);
    clear G R_eci R_ecef INTF;

    for s = 1:nS
        cfgS = cfg;
        cfgS.ffr.scheme = SCH{s}.scheme;
        cfgS.ffr.Delta  = SCH{s}.Delta;
        if ~isnan(SCH{s}.alpha), cfgS.ffr.alpha = SCH{s}.alpha; end

        A = ffr_allocate(cfgS, ctx, false);
        F = compute_sinr_ffr(cfgS, ctx, A);

        acc(s).SINR  = [acc(s).SINR,  F.SINR_dB];
        acc(s).CN    = [acc(s).CN,    F.CN_dB];
        acc(s).pen   = [acc(s).pen,   F.penalty_dB];
        acc(s).Bu    = [acc(s).Bu,    A.B_user_Hz];
        acc(s).isC   = [acc(s).isC,   A.isCenter];
        acc(s).alpha = [acc(s).alpha, A.alpha];
        % Interferencia INTER-satelite frente al ruido: es lo que decide si añadir
        % satelites cambia realmente algo (tarea 4). I_inter << N => no cambia.
        acc(s).Iint  = [acc(s).Iint,  F.I_inter_dBW];
        acc(s).Nn    = [acc(s).Nn,    F.N_dBW];

        % La SINR de referencia en reuso-1 es independiente del esquema: se cachea
        % dentro del bloque para que los dos esquemas clasifiquen igual.
        if isempty(ctx.SINR_ref_dB), ctx.SINR_ref_dB = A.SINR_ref_dB; end
    end

    elevAll = [elevAll, ctx.elev_deg];    %#ok<AGROW>
    covAll  = [covAll,  ctx.cov];         %#ok<AGROW>
    nVisAll = [nVisAll, ctx.nVis];        %#ok<AGROW>
    clear ctx;
end

E.secs         = toc(t0);
E.peak_abs_GB  = proc_peak_GB();                  % pico ABSOLUTO del proceso
E.peak_meas_GB = max(E.peak_abs_GB - pk0, 0);     % incremento atribuible al caso

% --- KPIs sobre las matrices concatenadas (identico a una pasada entera) ---
maskV = covAll & elevAll >= ST.elMin_valid;
E.K  = cell(1,nS);  E.KV = cell(1,nS);
for s = 1:nS
    FF = struct('SINR_dB',acc(s).SINR, 'CN_dB',acc(s).CN, 'penalty_dB',acc(s).pen);
    AL = struct('B_user_Hz',acc(s).Bu, 'isCenter',acc(s).isC, ...
                'scheme',SCH{s}.scheme, 'Delta',SCH{s}.Delta, 'alpha',acc(s).alpha);
    cfgS = cfg;  cfgS.ffr.scheme = SCH{s}.scheme;  cfgS.ffr.Delta = SCH{s}.Delta;
    E.K{s}  = strip_cdf(compute_kpis(cfgS, FF, AL));
    E.KV{s} = strip_cdf(compute_kpis(cfgS, FF, AL, maskV));

    % I_inter / N medio (en LINEAL, luego a dB): cuanto pesa la interferencia
    % INTER-satelite frente al ruido termico. Es la magnitud que decide si subir
    % la densidad de satelites cambia la interferencia o no.
    r = 10.^((acc(s).Iint - acc(s).Nn)/10);
    E.INR_dB(s) = 10*log10(mean(r(covAll & isfinite(r))));
end

E.covFrac   = sum(covAll(:)) / max(numel(covAll),1);
nv = nVisAll(logical(covAll));                 % logical(): defensa ante el tipo
E.nVis_mean = mean(nv(:));
end

% -------------------------------------------------------------------------
function gb = peak_model_GB(M, N, Nt, nB, ST)
%PEAK_MODEL_GB  Pico de memoria de UNA pasada entera (ver cabecera del script).
%   DELEGA en mem_peak_model_GB (FUENTE UNICA del modelo, compartida con
%   run_one_density y estimate_sweep_memory). Antes este script tenia su PROPIA
%   copia de la formula: es justo el patron "misma regla escrita en varios sitios"
%   que ya causo los bugs de merge_cfg/atm_key_of y del umbral de sesgo de layout.
%   Se conserva la firma con ST solo para no tocar los llamadores, y se COMPRUEBA
%   que las constantes del estudio siguen siendo las del modelo compartido.
if ST.bytes_per_elem ~= 74 || ST.bytes_per_belem ~= 40
    error('run_convergence_studyA:memModel', ...
        ['ST.bytes_per_elem/bytes_per_belem (%g/%g) ya no coinciden con ' ...
         'mem_peak_model_GB (74/40). Cambia el modelo en UN solo sitio.'], ...
        ST.bytes_per_elem, ST.bytes_per_belem);
end
gb = mem_peak_model_GB(M, N, Nt, nB);
end

% -------------------------------------------------------------------------
function [xc, dmax] = converged_step(xs, p5, tol, ok)
%CONVERGED_STEP  Valor de MENOR coste cuyo p5 ya coincide con TODOS los mas finos.
%   xs va ordenado de MENOS a MAS resolucion (step_km decreciente, o Nt creciente),
%   luego xs(i) "converge" si |p5(i) - p5(j)| < tol para TODO j > i. Comprobarlo
%   contra todos los refinamientos posteriores, y no solo contra el siguiente,
%   evita el falso positivo clasico: dos puntos consecutivos que coinciden por
%   casualidad en una curva que todavia no ha convergido.
xc = NaN;
n  = numel(xs);
dmax = nan(1, n-1);
for i = 1:n-1
    if ok(i) && ok(i+1)
        dmax(i) = max(abs(p5(:,i) - p5(:,i+1)));     % diagnostico: paso a paso
    end
end
for i = 1:n-1
    if ~ok(i), continue; end
    jj = find(ok((i+1):end)) + i;                    % todos los mas finos validos
    if isempty(jj), continue; end
    dAll = max(max(abs(p5(:,jj) - p5(:,i))));
    if dAll < tol
        xc = xs(i);
        return;
    end
end
end

% -------------------------------------------------------------------------
function gb = ram_free_GB()
gb = NaN;
try
    if ispc, [~, s] = memory; gb = s.PhysicalMemory.Available / 2^30; end
catch, gb = NaN;
end
if isnan(gb), gb = 4; end     % valor prudente si el sistema no lo expone
end

% -------------------------------------------------------------------------
function gb = ram_total_GB()
gb = NaN;
try
    if ispc, [~, s] = memory; gb = s.PhysicalMemory.Total / 2^30; end
catch, gb = NaN;
end
end

% -------------------------------------------------------------------------
function gb = proc_peak_GB()
%PROC_PEAK_GB  Pico de working set del proceso MATLAB [GB].
%   Es MONOTONO (maximo desde el arranque), asi que la diferencia antes/despues de
%   un caso solo es informativa si los casos se ejecutan de menor a mayor tamano,
%   que es como estan ordenados los barridos de este script.
try
    p = System.Diagnostics.Process.GetCurrentProcess();
    p.Refresh();
    gb = double(p.PeakWorkingSet64) / 2^30;
catch
    gb = NaN;
end
end

% -------------------------------------------------------------------------
function K = strip_cdf(K)
f = {'all','center','edge'};
for i = 1:numel(f)
    if isfield(K,f{i})
        K.(f{i}).cdf_SINR = [];  K.(f{i}).cdf_R = [];
    end
end
end

% -------------------------------------------------------------------------
function s = ternary_str(c, a, b)
if c, s = a; else, s = b; end
end

function s = ternax(c)
if c, s = 'FIABLE'; else, s = 'NO FIABLE'; end
end

% -------------------------------------------------------------------------
function [md, where] = deep_maxdiff(a, b, path)
%DEEP_MAXDIFF  Maxima diferencia absoluta recorriendo structs/cells (ver
%   test_parfor_invariance: misma logica, replicada aqui para que el script sea
%   autonomo).
md = 0;  where = path;
if ~strcmp(class(a), class(b))
    md = Inf;  where = sprintf('%s [clase]', path);  return;
end
if isstruct(a)
    if ~isequal(size(a), size(b)), md = Inf; where = [path ' [tam]']; return; end
    fa = sort(fieldnames(a));  fb = sort(fieldnames(b));
    if ~isequal(fa, fb), md = Inf; where = [path ' [campos]']; return; end
    for e = 1:numel(a)
        for i = 1:numel(fa)
            [d, w] = deep_maxdiff(a(e).(fa{i}), b(e).(fa{i}), sprintf('%s.%s', path, fa{i}));
            if d > md, md = d; where = w; end
        end
    end
    return;
end
if iscell(a)
    if ~isequal(size(a), size(b)), md = Inf; where = [path ' [tam cell]']; return; end
    for i = 1:numel(a)
        [d, w] = deep_maxdiff(a{i}, b{i}, sprintf('%s{%d}', path, i));
        if d > md, md = d; where = w; end
    end
    return;
end
if ischar(a) || isstring(a)
    if ~isequal(a,b), md = Inf; where = [path ' [texto]']; end
    return;
end
if isnumeric(a) || islogical(a)
    if ~isequal(size(a), size(b)), md = Inf; where = [path ' [tam num]']; return; end
    if isempty(a), return; end
    x = double(a);  y = double(b);
    d = abs(x-y);
    d(isnan(x) & isnan(y)) = 0;
    d(xor(isnan(x), isnan(y))) = Inf;
    md = max(d(:));
    return;
end
if ~isequaln(a,b), md = Inf; where = [path ' [otro]']; end
end
