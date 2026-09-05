%% RUN_CONVERGENCE_STEP_3SCH  Barrido de la REJILLA (step_km) con TRES esquemas.
%
%  HERMANO de run_convergence_nt_step4.m. No produce fisica nueva: da al barrido de
%  PASO DE REJILLA la misma cobertura de esquemas que ya tiene el barrido de Nt,
%  para que las dos tablas de convergencia se puedan presentar juntas.
%
%  EL HUECO. La tarea 1 de run_convergence_studyA (variable T1 de
%  convergence_studyA.mat) barre step_km = [8 6 4 3 2] con DOS esquemas
%  (reuse1 y ffr(Delta=3, alpha=0.4)). El barrido de Nt ya se ha ampliado a
%  ffr(Delta=4, alpha=0.4); aqui se hace lo mismo con el eje de rejilla.
%
%  QUE HACE
%    PARTE 1. Barrido de step_km con los TRES esquemas y, sobre esa MISMA tanda,
%      ANCLA contra T1 del .mat publicado: las filas 1:2 deben reproducirlo con
%      max|dif| = 0. Al comparar las dos filas publicadas DENTRO de un barrido de
%      tres, el ancla demuestra ademas que anadir ffr(D=4) no perturba a los
%      esquemas anteriores (los tres comparten el ctx dentro de eval_case).
%    PARTE 2. Aplica el criterio de convergencia YA EXISTENTE (sin tocarlo),
%      agregado sobre los 2 esquemas publicados (cifra historica, trazable) y
%      sobre los 3 actuales.
%    PARTE 3. Figura del barrido con las tres curvas.
%
%  ---------------------------------------------------------------------------
%  DEUDA TECNICA DECLARADA (leerla antes de tocar nada)
%  ---------------------------------------------------------------------------
%  eval_case, converged_step y sus auxiliares son COPIA LITERAL de las funciones
%  locales de run_convergence_studyA.m -- la MISMA copia que ya vive en
%  run_convergence_nt_step4.m, luego este fichero es la TERCERA. No es la fuente
%  unica que pide la politica del proyecto (una misma regla escrita en
%  varios sitios), pero las funciones LOCALES de un script no son invocables
%  desde fuera y extraerlas a ficheros propios obligaria a modificar
%  run_convergence_studyA.m -- un script ya publicado -- y a re-ejecutar la FASE A
%  entera, lo que sobrescribiria convergence_studyA.mat, que es justamente la
%  REFERENCIA del ancla.
%  MITIGACION: la PARTE 1 verifica la copia contra la SALIDA de la original. Si
%  alguien edita eval_case en run_convergence_studyA.m y no aqui, el ancla falla y
%  este script PARA. La divergencia silenciosa no es posible.
%
%  ---------------------------------------------------------------------------
%  MEMORIA
%  ---------------------------------------------------------------------------
%  El punto caro es step_km = 2 (M = 1257): el pico SIN TROCEAR modelado es
%  ~8.6 GB, mas que la RAM total de la maquina de referencia (7.45 GB). El
%  troceado temporal es EXACTO (test_timeblock_invariance / test_ctx_blocked:
%  max|dif| = 0) y lo acota a ~1.4 GB con bloques de 10 instantes.
%  OJO: el parametro de troceado de eval_case es su 4o argumento forceChunk, NO
%  cfg.compute.timeBlock. Se fija explicitamente (y eval_case toma el MINIMO entre
%  ese valor y lo que permita la RAM libre) para que el barrido sea reproducible.
%
%  NO toca el criterio (tol_dB = 0.2) ni ningun valor publicado.
%  Resultados en convergence_step_3sch.mat (NO sobrescribe convergence_studyA.mat).

clear; clc; close all;

%% 0. Parametros -- IDENTICOS a los de run_convergence_studyA (tarea 1)
ST.T_fixed   = 1584;
ST.P_fixed   = 72;
ST.radius_km = 40;
ST.window_s  = 7200;
ST.steps_km  = [8 6 4 3 2];          % mismo barrido que la tarea 1 original
ST.Nt_task1  = 61;                   % instantes FIJOS de la tarea 1 (dt = 120 s)
ST.tol_dB    = 0.2;                  % CRITERIO SIN TOCAR
ST.mem_frac        = 0.55;
ST.bytes_per_elem  = 74;
ST.bytes_per_belem = 40;
ST.elMin_valid     = 45;

FORCE_CHUNK  = 10;     % instantes por bloque  % <-- DECISION (ver cabecera)

%% 0b. Configuracion base -- COPIA LITERAL del cfg0 de run_convergence_studyA
cfg0 = config_default();
cfg0.ground.mode      = 'localgrid';
cfg0.ground.radius_km = ST.radius_km;
cfg0.ground.step_km   = 3;                 % lo sobrescribe cada caso
cfg0.time.t0          = 0;
cfg0.time.duration    = ST.window_s;
cfg0.ffr.tau_mode     = 'quantile';
cfg0.ffr.tau_q        = 50;
cfg0.ffr.sched        = 'share';
cfg0.beams.nRings     = 5;                 % 91 haces
cfg0.constellations(1).T = ST.T_fixed;
cfg0.constellations(1).P = ST.P_fixed;

%  Los DOS primeros esquemas son los publicados (y en ese orden): el ancla compara
%  las filas 1:2 contra T1, que se ejecuto con ellos dos.
%  'short' es SOLO una etiqueta de columna para las tablas estrechas.
SCH = { struct('scheme','reuse1','Delta',1,'alpha',NaN,'name','reuse1',        'short','r1'), ...
        struct('scheme','ffr',   'Delta',3,'alpha',0.4,'name','ffr(D=3,a=0.4)','short','ffr3'), ...
        struct('scheme','ffr',   'Delta',4,'alpha',0.4,'name','ffr(D=4,a=0.4)','short','ffr4') };
nSCH = numel(SCH);
nPUB = 2;                            % esquemas presentes en convergence_studyA.mat

% Contraste SECUNDARIO (informativo, no bloqueante) contra la tabla publicada del
% estudio de convergencia, en orden step = [8 6 4 3 2]. La puerta que MANDA es el
% ancla contra el .mat: una constante copiada a 3 decimales no se puede comprobar
% a 0 (su propio redondeo vale ya 5e-4). Misma politica que en h2_rerun.
PUB_TXT = [-4.997 -5.038 -4.950 -4.972 -4.974;    % reuse1
            3.645  3.757  3.865  3.782  3.812];   % ffr(D=3, a=0.4)
PUB_TXT_TOL = 6e-4;

fprintf('=========================================================================\n');
fprintf('  BARRIDO DE REJILLA (step_km) CON TRES ESQUEMAS\n');
fprintf('=========================================================================\n');
fprintf('T=%d (P=%d) | radio %d km | ventana %.0f h | Nt=%d (dt=%.0f s) | nRings=%d\n', ...
    ST.T_fixed, ST.P_fixed, ST.radius_km, ST.window_s/3600, ST.Nt_task1, ...
    ST.window_s/(ST.Nt_task1-1), cfg0.beams.nRings);
fprintf('step_km = %s\n', mat2str(ST.steps_km));
fprintf('Esquemas (%d):', nSCH);
for s = 1:nSCH, fprintf(' %s', SCH{s}.name); if s < nSCH, fprintf(' |'); end, end
fprintf('   (los %d primeros son los publicados)\n', nPUB);
fprintf('Troceado: forceChunk = %d instantes por bloque (EXACTO)\n', FORCE_CHUNK);
fprintf('RAM libre al arrancar: %.2f GB\n', ram_free_GB());

fprintf('\n[warmup] construyendo la tabla ITU-R P.618 (fuera del cronometro)...\n');
tw = tic;  atm_loss_dB(cfg0, 45);  tWarm = toc(tw);
fprintf('[warmup] listo (%.1f s)\n', tWarm);

tStudy = tic;

%% ------------------------------------------------------------------------
%  PARTE 1 -- Barrido de step_km (3 esquemas) + ANCLA contra T1 publicado
%  ------------------------------------------------------------------------
fprintf('\n=========== PARTE 1: barrido de step_km + ancla contra T1 ===========\n');
PUB = load('convergence_studyA.mat', 'T1');
assert(isequal(PUB.T1.step_km, ST.steps_km), ...
    'ANCLA: T1.step_km = %s no coincide con el barrido %s.', ...
    mat2str(PUB.T1.step_km), mat2str(ST.steps_km));
assert(size(PUB.T1.p5,1) == nPUB, ...
    'ANCLA: T1 del .mat trae %d esquemas, no %d.', size(PUB.T1.p5,1), nPUB);
assert(strcmp(SCH{1}.scheme,'reuse1') && SCH{1}.Delta == 1 && ...
       strcmp(SCH{2}.scheme,'ffr')    && SCH{2}.Delta == 3, ...
    ['ANCLA: los dos primeros esquemas de SCH ya no son (reuse1, ffr D=3), que ' ...
     'son los que ejecuto T1. El ancla compararia filas distintas.']);
fprintf('T1 publicado: step_km = %s, Nt = %d\n', mat2str(PUB.T1.step_km), ST.Nt_task1);
fprintf('NOTA: el original eligio nChunks = %s por RAM libre; aqui se fuerza %d.\n', ...
    mat2str(PUB.T1.nChunks), FORCE_CHUNK);
fprintf('      Si el troceado es exacto, el resultado debe ser IDENTICO igualmente.\n');
fprintf(['NOTA: el barrido corre con %d esquemas y el ancla compara SOLO las filas\n' ...
         '      1:%d (las publicadas). Que salga 0 prueba ademas que anadir\n' ...
         '      ffr(D=4) no perturba a los dos esquemas anteriores.\n\n'], nSCH, nPUB);

A = run_step_sweep(cfg0, SCH, ST, FORCE_CHUNK);

fprintf('\n%-8s %8s %10s', 'step_km','M','M_borde');
for s = 1:nSCH, fprintf(' %16s', SCH{s}.name); end
fprintf(' %8s %9s %8s\n', 'bloques','pico[GB]','t[s]');
fprintf('%s\n', repmat('-', 1, 28 + 17*nSCH + 28));
for i = 1:numel(ST.steps_km)
    fprintf('%-8g %8d %10d', ST.steps_km(i), A.M(i), A.Medge(i));
    for s = 1:nSCH, fprintf(' %+16.4f', A.p5(s,i)); end
    fprintf(' %8d %9.2f %8.1f\n', A.nChunks(i), A.peak_model_GB(i), A.secs(i));
end
fprintf('%s\n', repmat('-', 1, 28 + 17*nSCH + 28));
fprintf('(pico[GB] = pico modelado de una pasada ENTERA sin trocear; el troceado lo acota)\n');

% Tabla de GEOMETRIA VALIDA (elev >= 45 deg), informativa
fprintf('\n-- misma tabla restringida a geometria valida (elev >= %d deg) --\n', ST.elMin_valid);
fprintf('%-8s', 'step_km');
for s = 1:nSCH, fprintf(' %16s', SCH{s}.name); end
fprintf('\n%s\n', repmat('-', 1, 8 + 17*nSCH));
for i = 1:numel(ST.steps_km)
    fprintf('%-8g', ST.steps_km(i));
    for s = 1:nSCH, fprintf(' %+16.4f', A.p5V(s,i)); end
    fprintf('\n');
end

% --- ANCLA (bloqueante) ---
dAnc = max(max(abs(A.p5(1:nPUB,:) - PUB.T1.p5)));
fprintf('\nmax|dif| vs T1 publicado (filas 1:%d) = %.3e dB\n', nPUB, dAnc);
if dAnc ~= 0
    for s = 1:nPUB
        fprintf('  %-16s obtenido = %s\n', SCH{s}.name, mat2str(round(A.p5(s,:),6)));
        fprintf('  %-16s publicado= %s\n', '', mat2str(round(PUB.T1.p5(s,:),6)));
    end
    error('run_convergence_step_3sch:ancla', ...
        ['El ancla NO reproduce T1 (max|dif| = %.6g dB). O la copia de eval_case ' ...
         'ha divergido de la original, o el troceado no es exacto, o anadir el ' ...
         'tercer esquema ha perturbado a los anteriores. PARAR.'], dAnc);
end
fprintf('ANCLA OK: la copia de eval_case reproduce T1 EXACTAMENTE (max|dif| = 0).\n');

% --- Contraste secundario contra la tabla transcrita (informativo) ---
dTxt = max(max(abs(A.p5(1:nPUB,:) - PUB_TXT)));
fprintf('contraste vs tabla transcrita en CLAUDE.md (3 decimales): max|dif| = %.2e dB %s\n', ...
    dTxt, yesno_txt(dTxt <= PUB_TXT_TOL, '(dentro del redondeo)', ...
                    '(FUERA del redondeo: revisar la transcripcion)'));

%% ------------------------------------------------------------------------
%  PARTE 2 -- Criterio de convergencia (EL YA EXISTENTE, sin tocar)
%  ------------------------------------------------------------------------
fprintf('\n=========== PARTE 2: criterio de convergencia (tol = %.1f dB) ===========\n', ST.tol_dB);
fprintf('Regla (converged_step, copia literal): el step de MENOR coste (mas GRUESO)\n');
fprintf('cuyo p5 difiere menos de %.1f dB de TODAS las rejillas mas finas.\n\n', ST.tol_dB);

fprintf('%-8s', 'step_km');
for s = 1:nSCH, fprintf(' %14s', ['max|d| ' SCH{s}.short]); end
fprintf(' %10s %10s %8s\n', sprintf('max %d pub', nPUB), sprintf('max %d', nSCH), 'cumple');
fprintf('%s\n', repmat('-', 1, 8 + 15*nSCH + 32));
for i = 1:numel(ST.steps_km)
    fprintf('%-8g', ST.steps_km(i));
    if i < numel(ST.steps_km)
        d = nan(1,nSCH);
        for s = 1:nSCH
            d(s) = max(abs(A.p5(s,i+1:end) - A.p5(s,i)));
            fprintf(' %14.4f', d(s));
        end
        fprintf(' %10.4f %10.4f %8s\n', max(d(1:nPUB)), max(d), yesno(max(d) < ST.tol_dB));
    else
        for s = 1:nSCH, fprintf(' %14s', '-'); end
        fprintf(' %10s %10s %8s\n', '-','-','-');
    end
end

R = struct();
[R.step_conv_pub, R.dp5_pub] = converged_step(ST.steps_km, A.p5(1:nPUB,:), ST.tol_dB, A.ok);
[R.step_conv,     R.dp5    ] = converged_step(ST.steps_km, A.p5,           ST.tol_dB, A.ok);
fprintf('\ndifs consecutivas (max %d publicados): %s\n', nPUB, mat2str(round(R.dp5_pub,4)));
fprintf('difs consecutivas (max %d esquemas):   %s\n', nSCH, mat2str(round(R.dp5,4)));
fprintf('>> step_km CONVERGIDO con %d esquemas: %s\n', nPUB, conv_txt(R.step_conv_pub, ST.tol_dB));
fprintf('>> step_km CONVERGIDO con %d esquemas: %s\n', nSCH, conv_txt(R.step_conv, ST.tol_dB));

% La pregunta CONCRETA: el valor DE TRABAJO es step_km = 4. ¿Lo cumple cada esquema?
i4 = find(ST.steps_km == 4, 1);
R.d4_sch = max(abs(A.p5(:,i4+1:end) - A.p5(:,i4)), [], 2)';
R.d4     = max(R.d4_sch);
fprintf('\nPREGUNTA DIRECTA -- ¿el valor DE TRABAJO (step_km = 4) cumple el criterio?\n');
fprintf('  max|p5(step) - p5(4)| para rejillas mas finas, POR ESQUEMA:\n');
for s = 1:nSCH
    fprintf('    %-16s %.4f dB  ->  %s\n', SCH{s}.name, R.d4_sch(s), ...
        yesno_txt(R.d4_sch(s) < ST.tol_dB, 'CUMPLE', 'NO CUMPLE'));
end
fprintf('  agregado sobre los %d esquemas: %.4f dB (umbral %.1f dB)  ->  %s\n', ...
    nSCH, R.d4, ST.tol_dB, yesno_txt(R.d4 < ST.tol_dB, 'CUMPLE', 'NO CUMPLE'));
fprintf(['  RECORDATORIO: step_km = 4 es el valor DE TRABAJO por juicio de ingenieria\n' ...
         '  (mas fino que el convergido), no el resultado del criterio. Ver CLAUDE.md.\n']);

%% ------------------------------------------------------------------------
%  PARTE 3 -- Figura del barrido (las tres curvas)
%  ------------------------------------------------------------------------
%  DOS paneles porque la escala manda: los esquemas estan separados ~12 dB y en
%  valor absoluto la convergencia (decimas de dB) no se ve. El panel (b) da la
%  MAGNITUD QUE DECIDE, la del criterio: max|p5(step) - p5(step')| sobre TODAS las
%  rejillas mas finas.
FIGDIR = 'figs_faseA';
mk = {'o-','s-','^-'};
fF = figure('Name','Convergencia vs rejilla (3 esquemas)','Color','w','Visible','off', ...
            'Position',[100 100 1000 430]);

ax1 = subplot(1,2,1);  hold(ax1,'on');
for s = 1:nSCH
    plot(ax1, ST.steps_km(A.ok), A.p5(s,A.ok), mk{min(s,end)}, ...
        'LineWidth',1.6, 'DisplayName', SCH{s}.name);
end
xline(ax1, 4, 'k:', 'LineWidth',1.4, 'HandleVisibility','off');
set(ax1,'XDir','reverse'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'XTick', fliplr(ST.steps_km));
xlabel(ax1,'paso de la rejilla [km] (mas FINA hacia la derecha)');
ylabel(ax1,'SINR_{edge,p5} [dB]');
title(ax1, sprintf('(a) valor absoluto | T=%d, N_t=%d', ST.T_fixed, ST.Nt_task1));
legend(ax1,'Location','east');

% (b) la magnitud del criterio
dcrit = nan(nSCH, numel(ST.steps_km));
for i = 1:numel(ST.steps_km)-1
    dcrit(:,i) = max(abs(A.p5(:,i+1:end) - A.p5(:,i)), [], 2);
end
ax2 = subplot(1,2,2);  hold(ax2,'on');
for s = 1:nSCH
    plot(ax2, ST.steps_km(1:end-1), dcrit(s,1:end-1), mk{min(s,end)}, ...
        'LineWidth',1.6, 'DisplayName', SCH{s}.name);
end
yline(ax2, ST.tol_dB, 'k--', 'LineWidth',1.6, ...
      'Label',sprintf('criterio = %.1f dB', ST.tol_dB), 'FontSize',10, ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
xline(ax2, 4, 'k:', 'LineWidth',1.4, 'HandleVisibility','off');
set(ax2,'XDir','reverse'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'XTick', fliplr(ST.steps_km(1:end-1)));
ylim(ax2,[0, max(0.36, 1.15*max(dcrit(:)))]);
xlabel(ax2,'paso de la rejilla [km] (mas FINA hacia la derecha)');
ylabel(ax2,'max_{step'' mas fino} |\Delta p5| [dB]');
title(ax2,'(b) criterio de convergencia');
legend(ax2,'Location','northeast');

fpng = save_fig(fF, FIGDIR, 'faseA_e_convergencia_vs_paso_rejilla_3esquemas');
close(fF);
fprintf('\n[figura] %s\n', fpng);

%% 4. Guardar
timing = struct('warmup_p618_s', tWarm, 'study_s', toc(tStudy));
fprintf('\n[tiempos] barrido %.1f s (%.1f min) (+ %.1f s de warmup P.618)\n', ...
    timing.study_s, timing.study_s/60, timing.warmup_p618_s);
fprintf('[memoria] pico del proceso: %.2f GB | RAM libre ahora: %.2f GB\n', ...
    proc_peak_GB(), ram_free_GB());

save('convergence_step_3sch.mat', 'ST','cfg0','SCH','A','R', ...
     'FORCE_CHUNK','dAnc','dTxt','nSCH','nPUB','PUB_TXT','timing');
fprintf('\nResultados guardados en convergence_step_3sch.mat\n');
fprintf('(convergence_studyA.mat y convergence_nt_step4.mat NO se han tocado)\n');

% =========================================================================
% ============================ FUNCIONES ==================================
% =========================================================================

function A = run_step_sweep(cfg0, SCH, ST, forceChunk)
%RUN_STEP_SWEEP  Barrido de step_km con Nt fijo. Solo orquesta: la evaluacion la
%   hace eval_case, copia literal de run_convergence_studyA.
n  = numel(ST.steps_km);
nS = numel(SCH);
A = struct('step_km',ST.steps_km, 'Nt',ST.Nt_task1, 'M',nan(1,n), ...
           'Medge',nan(1,n), 'p5',nan(nS,n), 'p5V',nan(nS,n), ...
           'tens_GB',nan(1,n), 'peak_model_GB',nan(1,n), 'nChunks',nan(1,n), ...
           'ok',false(1,n), 'secs',nan(1,n));
for i = 1:n
    cfg = cfg0;
    cfg.ground.step_km = ST.steps_km(i);
    cfg.time.dt        = ST.window_s / (ST.Nt_task1 - 1);

    E = eval_case(cfg, SCH, ST, forceChunk);
    A.M(i)             = E.M;
    A.tens_GB(i)       = E.tens_GB;
    A.peak_model_GB(i) = E.peak_model_GB;
    A.nChunks(i)       = E.nChunks;
    A.secs(i)          = E.secs;
    A.ok(i)            = E.ok;
    if ~E.ok
        fprintf('  step=%-4g SALTADO: %s\n', ST.steps_km(i), E.why);
        continue;
    end
    A.Medge(i) = E.K{1}.viab.n_edge;
    for s = 1:nS
        A.p5(s,i)  = E.K{s}.viab.SINR_edge_p5;
        A.p5V(s,i) = E.KV{s}.viab.SINR_edge_p5;
    end
    fprintf('  step=%-4g M=%-5d ->', ST.steps_km(i), E.M);
    for s = 1:nS, fprintf(' %s %+8.4f |', SCH{s}.name, A.p5(s,i)); end
    fprintf('  (%d bloques, %.1f s)\n', E.nChunks, E.secs);
end
end

% -------------------------------------------------------------------------
function s = conv_txt(xc, tol)
if isnan(xc)
    s = sprintf('NO converge en el rango barrido con tol = %.1f dB', tol);
else
    s = sprintf('step_km = %g', xc);
end
end

% -------------------------------------------------------------------------
% A PARTIR DE AQUI: COPIA LITERAL de las funciones locales de
% run_convergence_studyA.m (las mismas que en run_convergence_nt_step4.m). No
% modificar sin modificar alli tambien; el ancla de la PARTE 1 detiene el script
% si divergen.
% -------------------------------------------------------------------------

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
        acc(s).Iint  = [acc(s).Iint,  F.I_inter_dBW];
        acc(s).Nn    = [acc(s).Nn,    F.N_dBW];

        if isempty(ctx.SINR_ref_dB), ctx.SINR_ref_dB = A.SINR_ref_dB; end
    end

    elevAll = [elevAll, ctx.elev_deg];    %#ok<AGROW>
    covAll  = [covAll,  ctx.cov];         %#ok<AGROW>
    nVisAll = [nVisAll, ctx.nVis];        %#ok<AGROW>
    clear ctx;
end

E.secs         = toc(t0);
E.peak_abs_GB  = proc_peak_GB();
E.peak_meas_GB = max(E.peak_abs_GB - pk0, 0);

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

    r = 10.^((acc(s).Iint - acc(s).Nn)/10);
    E.INR_dB(s) = 10*log10(mean(r(covAll & isfinite(r))));
end

E.covFrac   = sum(covAll(:)) / max(numel(covAll),1);
nv = nVisAll(logical(covAll));
E.nVis_mean = mean(nv(:));
end

% -------------------------------------------------------------------------
function gb = peak_model_GB(M, N, Nt, nB, ST)
if ST.bytes_per_elem ~= 74 || ST.bytes_per_belem ~= 40
    error('run_convergence_step_3sch:memModel', ...
        ['ST.bytes_per_elem/bytes_per_belem (%g/%g) ya no coinciden con ' ...
         'mem_peak_model_GB (74/40). Cambia el modelo en UN solo sitio.'], ...
        ST.bytes_per_elem, ST.bytes_per_belem);
end
gb = mem_peak_model_GB(M, N, Nt, nB);
end

% -------------------------------------------------------------------------
function [xc, dmax] = converged_step(xs, p5, tol, ok)
%CONVERGED_STEP  Valor de MENOR coste cuyo p5 ya coincide con TODOS los mas finos.
xc = NaN;
n  = numel(xs);
dmax = nan(1, n-1);
for i = 1:n-1
    if ok(i) && ok(i+1)
        dmax(i) = max(abs(p5(:,i) - p5(:,i+1)));
    end
end
for i = 1:n-1
    if ~ok(i), continue; end
    jj = find(ok((i+1):end)) + i;
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
if isnan(gb), gb = 4; end
end

% -------------------------------------------------------------------------
function gb = proc_peak_GB()
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
function s = yesno(b)
if b, s = 'SI'; else, s = 'NO'; end
end

function s = yesno_txt(c, a, b)
if c, s = a; else, s = b; end
end
