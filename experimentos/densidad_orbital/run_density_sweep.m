%% RUN_DENSITY_SWEEP  Eje de DENSIDAD ORBITAL del capitulo 6 (Sec. 6.2).
%  RUNNER DELGADO: no contiene fisica. Construye el perfil de estudio Ku, la lista
%  de puntos del barrido (densidad x esquema) y delega en el motor:
%
%      run_sweep_points  ->  parfor  ->  run_one_density  ->  pipeline completo
%
%  Cada punto de densidad es INDEPENDIENTE (run_one_density es una funcion pura),
%  luego el barrido se paraleliza sin reducciones ni estado compartido. La
%  equivalencia numerica con el barrido secuencial esta verificada en
%  test_parfor_invariance.
%
%  ---------------------------------------------------------------------------
%  QUE ES ESTE RUNNER Y QUE NO ES (leerlo antes de citar sus cifras)
%  ---------------------------------------------------------------------------
%  Este es el eje de densidad orbital DEFINITIVO del capitulo 6: barre la
%  constelacion completa (T y P a la vez, capa Walker realista) sobre EXACTAMENTE
%  el escenario de cabecera de run_ffr_demo, y con los CINCO esquemas de referencia
%  del capitulo mas la FFR adaptativa.
%
%  NO sustituye a run_convergence_studyA (FASE A), que barre el mismo eje T pero
%  con OTRO proposito: alli es un ESTUDIO DE CONVERGENCIA Y MEMORIA (rejilla de
%  6 km, Nt=61, dos esquemas, sin veredictos) que sirve para DIMENSIONAR los
%  experimentos. FASE A se queda como esta: es el precedente de dimensionado, y
%  este runner es el eje del capitulo.
%
%  SALIDA DE DECISION (criterio de viabilidad, definido en compute_kpis):
%  KPI.viab.verdict en cada punto y para cada esquema. La densidad MAXIMA con
%  verdict = 'viable' es el MARGEN DE VIABILIDAD, y su desplazamiento de reuso-1 a
%  FFR es el resultado central de H1. Se registra tambien margin_floor_dB (distancia
%  continua al suelo fisico, permite interpolar la densidad de cruce entre puntos)
%  y el criterio de H2  Delta*SE_centro/SE_borde  en cada densidad: si al saturar
%  cruza 1, la FFR adaptativa recuperaria margen y el resultado negativo de E5
%  seria especifico de la densidad baja (T=66), no general.
%
%  LAS DOS LECTURAS (K y KV) SE CONSERVAN A PROPOSITO. K son los KPIs de la
%  VENTANA COMPLETA y KV los de la submuestra de GEOMETRIA VALIDA (elev >= 45 deg).
%  Su CONTRASTE es el resultado, no un detalle de presentacion: si el margen de
%  viabilidad cambia mucho entre las dos columnas, la frontera la esta fijando la
%  ELEVACION del servidor y no la densidad de satelites.
%
%  ESQUEMAS COMPARADOS: los CINCO del capitulo 6 (los mismos, con los mismos Delta
%  y el mismo alpha = 0.4, que run_ffr_demo) MAS la FFR ADAPTATIVA como sexto caso.
%
%  Uso:  run_density_sweep        (ajustar antes Tlist / nWorkers segun la maquina)
%  Resultados en e3_density_sweep.mat
%
%  AVISO DE COSTE: cada worker mantiene SU PROPIA geometria (M*N*Ntb). Ver el
%  bloque de troceado mas abajo: run_sweep_points imprime la estimacion antes de
%  arrancar; si no cabe, bajar nWorkers, la rejilla o timeBlock ANTES de lanzar.

clear; clc; close all;

%% 1. Perfil de estudio (identico a run_ffr_demo salvo la densidad, que se barre)
cfg = config_default();

% REJILLA Y VENTANA: son DELIBERADAMENTE las de run_ffr_demo (E2), no un ajuste de
% coste de este runner. Rejilla local de radio 40 km y paso 3 km -> M = 553
% usuarios; ventana de 2 h con dt = 60 s -> Nt = 121 instantes. Se eligen asi para
% que el punto T = 1584 / P = 72 de este barrido sea EL MISMO PUNTO que E2 y deba
% por tanto REPRODUCIR ffr_results.mat exactamente. Esa reproduccion se comprueba
% en la Sec. 9 y es lo que ancla todo el eje: sin ella, el barrido daria una curva
% plausible pero sin punto de contraste con el escenario de cabecera del capitulo.
% (Ojo: la escalera de pasos del proyecto es 3 km en E2/E5, 4 km en los barridos y
% 6 km el convergido de FASE A; aqui se usa 3 km POR EL ANCLA, no por convergencia.)
cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;           % = run_ffr_demo
cfg.ground.step_km   = 3;            % = run_ffr_demo  (M = 553)

cfg.time.dt       = 60;              % s   = run_ffr_demo  (Nt = 121)
cfg.time.duration = 7200;            % s (2 h) = run_ffr_demo

cfg.ffr.tau_mode = 'quantile';   % reparto fijo 50/50 (ver run_ffr_demo)
cfg.ffr.tau_q    = 50;
cfg.ffr.sched    = 'share';

cfg.beams.nRings = 5;            % 91 haces (>= 1.5 anillos de guarda co-canal)

% TROCEADO TEMPORAL: OBLIGATORIO en la maquina de referencia, y el valor lo fija
% el PEOR punto del eje, que ahora es T = 4000 (no T = 1584). Con el modelo de
% memoria del proyecto (mem_peak_model_GB, cota superior de la FASE A) el pico de
% UN bloque es
%     M*Ntb*(N*74 + nBeams*40)/2^30
% y con M=553, nBeams=91:
%     T=1584, Ntb=10 -> 0.62 GB   (es el pico medido de ffr_results.mat: 0.6224 GB)
%     T=4000, Ntb=10 -> 1.54 GB   <- POR ENCIMA de la RAM libre tipica (~1.9 GB)
%     T=4000, Ntb= 4 -> 0.62 GB   <- mismo pico que el punto ya ejecutado sin
%                                    incidencias en E2, luego es un valor seguro
% Se fija por tanto timeBlock = 4. El troceado es EXACTO (el pipeline no acopla
% instantes; test_timeblock_invariance da max|dif| = 0), luego cambiar 10 -> 4 no
% altera ningun numero: solo acota el pico. De hecho la comprobacion de la Sec. 9
% contrasta este barrido (Ntb=4) contra ffr_results.mat (Ntb=10) y debe dar 0.
cfg.compute.timeBlock = 4;                                         % <-- DECISION

elMin_valid = 45;                % deg  submuestra de GEOMETRIA VALIDA

%% 2. Rejilla de densidades: se escalan T Y P (capa Walker realista)  % <-- DECISION
%  Se barre la constelacion COMPLETA, escalando el numero de planos P junto con el
%  numero total de satelites T de modo que la razon T/P (satelites por plano) se
%  mantenga entre 22 y 25, que es el regimen de las capas Walker reales
%  (Starlink Shell-1 es 1584/72 = 22).
%
%  POR QUE NO SE DEJA P FIJO: la version anterior de este runner barria T con
%  P = 6 planos fijos. Eso hace que el punto T = 1584 sean 264 satelites POR PLANO
%  -- seis planos densisimos en vez de una capa Walker -- con el resultado de que
%  el servidor queda sistematicamente a elevacion baja (elev media 51.1 deg,
%  minima 25.7) y reuso-1 sale 'inviable' a -11.07 dB, frente a los -5.04 dB que
%  da la MISMA densidad con P = 72. La diferencia no es fisica: es la geometria de
%  la capa. Las listas de abajo son las de la FASE A, que ya escalaba P.
Tlist = [66 300 600 1000 1584 2500 4000];
Plist = [ 3  12  24   40   72  100  160];

% Guarda: una constelacion Walker exige T multiplo de P (T/P satelites por plano).
if numel(Tlist) ~= numel(Plist)
    error('run_density_sweep:listas', ...
        'Tlist (%d) y Plist (%d) deben tener la misma longitud.', numel(Tlist), numel(Plist));
end
bad = mod(Tlist, Plist) ~= 0;
if any(bad)
    error('run_density_sweep:walker', ...
        ['T debe ser multiplo de P en una constelacion Walker. Fallan: %s ' ...
         '(T/P = %s).'], mat2str(Tlist(bad)), mat2str(Tlist(bad)./Plist(bad)));
end
fprintf('Eje de densidad orbital (Walker, P escalado con T):\n');
for i = 1:numel(Tlist)
    fprintf('   T = %5d | P = %4d | %5.1f satelites/plano\n', ...
        Tlist(i), Plist(i), Tlist(i)/Plist(i));
end

%% 3. Esquemas por punto: los CINCO del capitulo 6 + la FFR ADAPTATIVA
%  Los cinco primeros son EXACTAMENTE los de run_ffr_demo (mismo scheme, mismo
%  Delta, mismo alpha y el MISMO NOMBRE), porque el punto T=1584/P=72 debe
%  reproducirlos uno a uno en la Sec. 9. El sexto (adaptativa) no esta en E2 y no
%  entra en esa comprobacion, pero se evalua aqui porque su crit_ratio_ffr a lo
%  largo del eje de densidad es material del capitulo 7 (H2) y no puede perderse:
%  es lo que demuestra que la refutacion de H2 no es un artefacto de T = 66.
%  Todos comparten geometria, instantes y clasificacion centro/borde de referencia
%  (run_one_density aborta si esta ultima difiere entre esquemas).
ALPHA = 0.4;                     % = run_ffr_demo (NO el alpha* = 0.5 de E3/E3b/E4)
schemeCases = { ...
    struct('scheme','reuse1','Delta',1,'alpha',NaN,  'adaptive',false, 'name','reuse1'), ...
    struct('scheme','reuseD','Delta',3,'alpha',NaN,  'adaptive',false, 'name','reuseD (D=3)'), ...
    struct('scheme','reuseD','Delta',4,'alpha',NaN,  'adaptive',false, 'name','reuseD (D=4)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',ALPHA,'adaptive',false, 'name','FFR (D=3, a=0.4)'), ...
    struct('scheme','ffr',   'Delta',4,'alpha',ALPHA,'adaptive',false, 'name','FFR (D=4, a=0.4)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',ALPHA,'adaptive',true,  'name','FFR-ADAPTATIVA (D=3)') };
nS = numel(schemeCases);
nRef = 5;                        % los 5 primeros son los comunes con run_ffr_demo

%% 4. Lista de puntos
pts = cell(1, numel(Tlist));
for i = 1:numel(Tlist)
    pts{i} = struct( ...
        'label',       sprintf('T=%d', Tlist(i)), ...
        'T',           Tlist(i), ...
        'P',           Plist(i), ...
        'cases',       {schemeCases}, ...
        'elMin_valid', elMin_valid, ...
        'keep_cdf',    false, ...   % las CDF se reconstruyen reejecutando un punto
        'seed',        42);
end

%% 5. Barrido. nWorkers = [] -> todos los disponibles; 0 -> secuencial.
%  Se fija 0 (SECUENCIAL) en la maquina de referencia: aun con troceado cada worker
%  necesita su propia geometria (0.62 GB en el peor punto) y en un portatil de
%  7.4 GB no caben los 4. No se pierde nada metodologicamente: parfor con
%  nWorkers=0 corre POR EL MISMO CAMINO DE CODIGO, y la equivalencia numerica
%  secuencial/paralelo esta verificada aparte en test_parfor_invariance
%  (isequaln = si, max|dif| = 0 exacto). En una maquina con mas RAM, poner [].
nWorkers = 0;                                                      % <-- DECISION
[results, INFO] = run_sweep_points(cfg, pts, nWorkers);

%% 6. Tablas del CRITERIO DE VIABILIDAD, en las DOS lecturas
%  Se imprimen las dos a proposito (ver cabecera): el contraste entre la ventana
%  completa y la geometria valida es lo que revela si la frontera del margen la fija
%  la ELEVACION o la densidad.
print_verdict_table(results, schemeCases, 'K',  'VENTANA COMPLETA', elMin_valid);
print_verdict_table(results, schemeCases, 'KV', 'GEOMETRIA VALIDA', elMin_valid);

%% 7. MARGEN DE VIABILIDAD: densidad maxima con verdict = 'viable' por esquema
[Tmax,  ok_K]  = margin_of(results, nS, 'K');
[TmaxV, ok_KV] = margin_of(results, nS, 'KV');

fprintf('\n[H1] MARGEN DE VIABILIDAD (T maxima con veredicto ''viable''):\n');
fprintf('   %-26s %-28s %-28s\n', 'esquema', 'ventana completa', sprintf('geometria valida (>=%d deg)', elMin_valid));
fprintf('   %s\n', repmat('-', 1, 84));
for c = 1:nS
    fprintf('   %-26s %-28s %-28s\n', schemeCases{c}.name, ...
        margin_str(Tmax(c)), margin_str(TmaxV(c)));
end
fprintf(['\n   Lectura: si un esquema es ''viable'' en geometria valida pero no en la\n' ...
         '   ventana completa, lo que le falta NO es densidad de satelites sino\n' ...
         '   ELEVACION del servidor: la frontera del margen la fija la geometria.\n']);

%% 8. Criterio de H2 por densidad:  Delta*SE_centro/SE_borde > 1
%  El motor lo devuelve NaN fuera de los esquemas 'ffr' (unicos con sub-banda
%  INTERIOR); en reuse1/reuseD el cociente no significa nada. Se imprime "-".
fprintf('\n[H2] Criterio Delta*SE_centro/SE_borde (la sub-banda interior renta si > 1):\n');
fprintf('     "-" = el criterio NO aplica a ese esquema (no tiene sub-banda interior).\n');
fprintf('%-8s', 'T');
for c = 1:nS, fprintf(' %22s', schemeCases{c}.name); end
fprintf('\n');
for i = 1:numel(results)
    fprintf('%-8d', results{i}.T);
    for c = 1:nS
        cr = results{i}.diag.crit_ratio_ffr(c);
        if isnan(cr), fprintf(' %22s', '-'); else, fprintf(' %22.3f', cr); end
    end
    fprintf('\n');
end

%% 8b. Diagnostico geometrico del eje (lo que hace que el eje sea el que es)
fprintf('\n[geometria] cobertura y servidor por densidad:\n');
fprintf('%-8s %6s %8s %10s %10s %10s %10s\n', 'T','P','covFrac','elev_min','elev_p5','elev_mean','nVis');
for i = 1:numel(results)
    d = results{i}.diag;
    fprintf('%-8d %6d %8.4f %10.2f %10.2f %10.2f %10.3f\n', ...
        results{i}.T, results{i}.P, d.covFrac, d.elev_min, d.elev_p5, d.elev_mean, d.nVis_mean);
end

%% 9. VERIFICACION OBLIGATORIA: el punto T=1584/P=72 debe reproducir ffr_results.mat
%  Es el ancla del eje. Este barrido llega al mismo punto por OTRO camino de codigo
%  (run_one_density + troceado con timeBlock=4 + bucle de barrido) que el de E2
%  (build_ctx_blocked con timeBlock=10), asi que la coincidencia valida a la vez la
%  funcion de punto, el troceado y el bucle. Se compara la VENTANA COMPLETA (K) del
%  punto contra KK de ffr_results.mat en los CINCO esquemas comunes.
fprintf('\n========== VERIFICACION: T=1584/P=72 vs ffr_results.mat (E2) ==========\n');
iRef = find(Tlist == 1584 & Plist == 72, 1);
maxdif = NaN;  nMismatch = 0;  okRef = false;
if isempty(iRef)
    warning('run_density_sweep:sinAncla', ...
        'El eje no contiene el punto T=1584/P=72: NO se puede verificar contra E2.');
elseif ~isfile('ffr_results.mat')
    warning('run_density_sweep:sinE2', ...
        'No se encuentra ffr_results.mat: NO se puede verificar el ancla. Ejecutar run_ffr_demo.');
else
    E2 = load('ffr_results.mat', 'cases', 'KK');
    r  = results{iRef};
    fprintf('%-20s %14s %14s %14s | %s\n', 'campo', 'barrido', 'E2', '|dif|', 'esquema');
    maxdif = 0;
    for c = 1:nRef
        % Emparejamiento por NOMBRE (no por indice): si alguien reordena los casos
        % de un runner, esto lo detecta en vez de comparar esquemas distintos.
        j = find(strcmp(cellfun(@(x) x.name, E2.cases, 'uni', 0), schemeCases{c}.name), 1);
        if isempty(j)
            warning('run_density_sweep:sinEsquema', ...
                'El esquema "%s" no existe en ffr_results.mat: no se compara.', schemeCases{c}.name);
            nMismatch = nMismatch + 1;
            continue;
        end
        A = r.K{c};   B = E2.KK{j};
        vals = { 'SINR_edge_p5', A.viab.SINR_edge_p5,   B.viab.SINR_edge_p5;   ...
                 'R_p5_borde',   A.edge.R_p5_Mbps,      B.edge.R_p5_Mbps;      ...
                 'R_agg',        A.all.R_agg_Mbps,      B.all.R_agg_Mbps;      ...
                 'n_edge',       A.viab.n_edge,         B.viab.n_edge };
        for k = 1:size(vals,1)
            d = abs(vals{k,2} - vals{k,3});
            maxdif = max(maxdif, d);
            fprintf('%-20s %14.6f %14.6f %14.3e | %s\n', ...
                vals{k,1}, vals{k,2}, vals{k,3}, d, schemeCases{c}.name);
        end
        if ~strcmpi(A.viab.verdict, B.viab.verdict)
            nMismatch = nMismatch + 1;
            fprintf('  !! VEREDICTO DISTINTO en %s: barrido="%s" vs E2="%s"\n', ...
                schemeCases{c}.name, A.viab.verdict, B.viab.verdict);
        end
    end
    okRef = (maxdif == 0) && (nMismatch == 0);
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('max|dif| = %.3e | veredictos discrepantes: %d | %s\n', maxdif, nMismatch, ...
        ternary_str(okRef, 'ANCLA REPRODUCIDA (exacto)', 'DISCREPANCIA'));
end

if ~okRef
    % No se sobreescribe el .mat publicado: el resultado se aparta en cuarentena
    % para poder diagnosticar sin perder el barrido, y se aborta.
    save('e3_density_sweep_FAILED.mat', 'cfg', 'pts', 'results', 'INFO', ...
         'Tlist', 'Plist', 'schemeCases', 'elMin_valid', '-v7.3');
    error('run_density_sweep:ancla', ...
        ['El punto T=1584/P=72 NO reproduce ffr_results.mat (max|dif| = %.3e, ' ...
         '%d veredictos discrepantes). NO se ha guardado e3_density_sweep.mat; el ' ...
         'barrido queda en e3_density_sweep_FAILED.mat para diagnostico. Causas a ' ...
         'descartar por orden: (1) ffr_results.mat corresponde a otra configuracion ' ...
         '(comprobar step_km=3, dt=60, nRings=5, tau_q=50, alpha=0.4); (2) alguien ' ...
         'ha tocado el motor desde que se genero E2; (3) el troceado ha dejado de ' ...
         'ser exacto (relanzar test_timeblock_invariance).'], maxdif, nMismatch);
end

%% 10. Guardar
save('e3_density_sweep.mat', 'cfg', 'pts', 'results', 'INFO', 'Tlist', 'Plist', ...
     'schemeCases', 'elMin_valid', 'Tmax', 'TmaxV', 'ok_K', 'ok_KV', ...
     'maxdif', 'nRef', '-v7.3');
fprintf('\nResultados guardados en e3_density_sweep.mat (%.1f s, %s)\n', ...
    INFO.elapsed_s, INFO.mode);

% =========================================================================
function print_verdict_table(results, schemeCases, campo, titulo, elMin_valid)
%PRINT_VERDICT_TABLE  Tabla SINR_edge_p5 + veredicto por densidad y esquema.
nS = numel(schemeCases);
if strcmp(campo,'KV')
    fprintf('\n========== VEREDICTO DE VIABILIDAD | %s (elev >= %d deg) ==========\n', ...
        titulo, elMin_valid);
else
    fprintf('\n========== VEREDICTO DE VIABILIDAD | %s ==========\n', titulo);
end
fprintf('%-8s %6s %6s', 'T', 'P', 'N');
for c = 1:nS, fprintf(' %24s', schemeCases{c}.name); end
fprintf('\n%s\n', repmat('-', 1, 22 + 25*nS));
for i = 1:numel(results)
    r = results{i};
    K = r.(campo);
    fprintf('%-8d %6d %6d', r.T, r.P, r.N);
    for c = 1:nS
        if isempty(K) || isempty(K{c})
            fprintf(' %24s', '(sin datos)');
        else
            V = K{c}.viab;
            fprintf(' %11.2f dB %-10s', V.SINR_edge_p5, upper(V.verdict));
        end
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('-', 1, 22 + 25*nS));
end

% -------------------------------------------------------------------------
function [Tmax, ok] = margin_of(results, nS, campo)
%MARGIN_OF  Densidad maxima con veredicto 'viable' por esquema.
Tmax = nan(1, nS);
ok   = false(nS, numel(results));
for c = 1:nS
    for i = 1:numel(results)
        K = results{i}.(campo);
        if isempty(K) || isempty(K{c}), continue; end
        ok(c,i) = strcmpi(K{c}.viab.verdict, 'viable');
    end
    if any(ok(c,:))
        Ts = cellfun(@(r) r.T, results);
        Tmax(c) = max(Ts(ok(c,:)));
    end
end
end

% -------------------------------------------------------------------------
function s = margin_str(T)
if isnan(T), s = 'ninguna densidad viable'; else, s = sprintf('%d satelites', T); end
end

% -------------------------------------------------------------------------
function s = ternary_str(c, a, b)
if c, s = a; else, s = b; end
end
