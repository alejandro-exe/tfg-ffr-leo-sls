%% RUN_CONVERGENCE_NT_STEP4  Barrido de Nt con la rejilla de TRABAJO (step_km = 4).
%
%  CIERRA UN HUECO DE VALIDACION DE LA FASE A, no produce fisica nueva.
%
%  EL HUECO. La tarea 2 de run_convergence_studyA (variable T2 de
%  convergence_studyA.mat) barre Nt con el step CONVERGIDO de la tarea 1, porque
%  run_convergence_studyA.m:215 hace  step_use = step_conv  y  step_conv = 6.
%  Pero el valor de TRABAJO de todos los experimentos del proyecto es
%  step_km = 4 (M = 317, no M = 137). Es decir: Nt = 61 se valido bajo una
%  rejilla distinta de aquella con la que despues se usa. Es la limitacion
%  estandar de un estudio de convergencia ENCADENADO (cada eje se barre fijando
%  el anterior en su valor convergido); aqui se MIDE en lugar de argumentarse.
%
%  QUE HACE
%    PARTE 1 (ANCLA). Repite el barrido a step_km = 6 y comprueba que reproduce
%      T2 del .mat publicado con max|dif| = 0. Es lo que legitima la comparacion:
%      sin ella, cualquier diferencia de la PARTE 2 seria indistinguible de un
%      artefacto de haber copiado mal la funcion de evaluacion.
%    PARTE 2. Mismo barrido con step_km = 4, TODO lo demas identico.
%    PARTE 3. Aplica el criterio de convergencia YA EXISTENTE (sin tocarlo) y
%      compara las dos rejillas punto a punto.
%
%  ---------------------------------------------------------------------------
%  AMPLIACION A TRES ESQUEMAS (peticion del tutor)
%  ---------------------------------------------------------------------------
%  El barrido incluye ahora ffr(Delta = 4, alpha = 0.4) ademas de reuse1 y
%  ffr(Delta = 3, alpha = 0.4). Es una AMPLIACION, no un cambio: eval_case ya era
%  generico en numel(SCH) y no se ha tocado; lo unico que estaba escrito a dos
%  esquemas era el orquestador (run_nt_sweep) y las tablas de impresion.
%    - El ANCLA (PARTE 1) sigue comparando SOLO las dos filas publicadas
%      (reuse1 y ffr(D=3)) contra T2, que se ejecuto con dos esquemas. Al hacerlo
%      sobre un barrido de TRES, el ancla verifica ademas que anadir el tercer
%      esquema NO perturba a los dos anteriores (comparten ctx dentro de
%      eval_case), que es justo lo que habria que demostrar.
%    - La PARTE 5 vuelve a comprobar lo mismo en la rejilla de TRABAJO, contra los
%      valores publicados del estudio. Si alguno se moviese, el script PARA.
%  El criterio (tol_dB = 0.2) NO se toca; se reporta agregado sobre los DOS
%  esquemas publicados y sobre los TRES, para que la cifra antigua siga siendo
%  trazable.
%
%  ---------------------------------------------------------------------------
%  DEUDA TECNICA DECLARADA (leerla antes de tocar nada)
%  ---------------------------------------------------------------------------
%  eval_case, converged_step y sus auxiliares son COPIA LITERAL de las funciones
%  locales de run_convergence_studyA.m. No es la fuente unica que pide la politica
%  del proyecto (una misma regla escrita en varios sitios), pero las
%  funciones LOCALES de un script no son invocables desde fuera y extraerlas a
%  ficheros propios obligaria a modificar un script ya publicado y a re-ejecutar
%  FASE A entera para demostrar que nada cambia.
%  MITIGACION: la PARTE 1 verifica la copia contra la SALIDA de la original. Si
%  alguien edita eval_case en run_convergence_studyA.m y no aqui, el ancla falla
%  y este script PARA. La divergencia silenciosa no es posible.
%
%  ---------------------------------------------------------------------------
%  MEMORIA
%  ---------------------------------------------------------------------------
%  A step_km = 4 (M = 317) con Nt = 241 el pico SIN TROCEAR es ~8.6 GB, mas que
%  la RAM total de la maquina de referencia (7.45 GB). El troceado temporal es
%  EXACTO (test_timeblock_invariance / test_ctx_blocked: max|dif| = 0).
%  OJO: el parametro de troceado de eval_case es su 4o argumento forceChunk, NO
%  cfg.compute.timeBlock (ese campo lo lee run_one_density, que es OTRO camino de
%  codigo y no interviene aqui). Se fija explicitamente en vez de dejar la
%  eleccion automatica por RAM libre de eval_case, para que el barrido sea
%  REPRODUCIBLE y no dependa de cuanta memoria hubiera libre ese dia.
%
%  NO toca el criterio (tol_dB = 0.2) ni ningun valor publicado.
%  Resultados en convergence_nt_step4.mat (NO sobrescribe convergence_studyA.mat).

clear; clc; close all;

%% 0. Parametros -- IDENTICOS a los de run_convergence_studyA
ST.T_fixed   = 1584;
ST.P_fixed   = 72;
ST.radius_km = 40;
ST.window_s  = 7200;
ST.Nt_list   = [31 61 121 241];      % mismo barrido que la tarea 2 original
ST.tol_dB    = 0.2;                  % CRITERIO SIN TOCAR
ST.mem_frac        = 0.55;
ST.bytes_per_elem  = 74;
ST.bytes_per_belem = 40;
ST.elMin_valid     = 45;

STEP_ANCLA   = 6;      % la rejilla con la que se ejecuto T2 (para el ancla)
STEP_TRABAJO = 4;      % la rejilla de trabajo de E2/E3/E3b/E4/EB/nRings
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

%  Los DOS primeros esquemas son los publicados (y en ese orden): el ancla de la
%  PARTE 1 compara las filas 1:2 contra T2, que se ejecuto con ellos dos.
%  'short' es SOLO una etiqueta de columna para las tablas estrechas.
SCH = { struct('scheme','reuse1','Delta',1,'alpha',NaN,'name','reuse1',        'short','r1'), ...
        struct('scheme','ffr',   'Delta',3,'alpha',0.4,'name','ffr(D=3,a=0.4)','short','ffr3'), ...
        struct('scheme','ffr',   'Delta',4,'alpha',0.4,'name','ffr(D=4,a=0.4)','short','ffr4') };
nSCH   = numel(SCH);
nPUB   = 2;                          % esquemas presentes en convergence_studyA.mat

% Valores PUBLICADOS a step_km = 4 (Nt = 61 re-verificado con la rejilla de
% TRABAJO), en orden Nt = [31 61 121 241]. Estan dados a 4 decimales,
% luego la tolerancia de la comprobacion es el propio redondeo (5e-5).
PUB_STEP4 = [-5.0080 -4.9504 -4.9461 -4.9066;    % reuse1
              3.6635  3.8653  3.8788  3.9739];   % ffr(D=3, a=0.4)
PUB_TOL   = 6e-5;

fprintf('=========================================================================\n');
fprintf('  BARRIDO DE Nt CON LA REJILLA DE TRABAJO (step_km = %d)\n', STEP_TRABAJO);
fprintf('=========================================================================\n');
fprintf('T=%d (P=%d) | radio %d km | ventana %.0f h | nRings=%d | tau cuantil %d%%\n', ...
    ST.T_fixed, ST.P_fixed, ST.radius_km, ST.window_s/3600, ...
    cfg0.beams.nRings, cfg0.ffr.tau_q);
fprintf('Nt = %s  ->  dt = %s s\n', mat2str(ST.Nt_list), ...
    mat2str(ST.window_s./(ST.Nt_list-1)));
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
%  PARTE 1 -- ANCLA: reproducir T2 (step_km = 6) del .mat publicado
%  ------------------------------------------------------------------------
fprintf('\n=========== PARTE 1: ANCLA contra convergence_studyA.mat (T2) ===========\n');
PUB = load('convergence_studyA.mat', 'T2');
assert(PUB.T2.step_km == STEP_ANCLA, ...
    'ANCLA: T2 del .mat se ejecuto con step_km = %g, no %g.', PUB.T2.step_km, STEP_ANCLA);
assert(isequal(PUB.T2.Nt, ST.Nt_list), ...
    'ANCLA: T2.Nt = %s no coincide con el barrido %s.', ...
    mat2str(PUB.T2.Nt), mat2str(ST.Nt_list));
assert(size(PUB.T2.p5,1) == nPUB, ...
    'ANCLA: T2 del .mat trae %d esquemas, no %d.', size(PUB.T2.p5,1), nPUB);
assert(strcmp(SCH{1}.scheme,'reuse1') && SCH{1}.Delta == 1 && ...
       strcmp(SCH{2}.scheme,'ffr')    && SCH{2}.Delta == 3, ...
    ['ANCLA: los dos primeros esquemas de SCH ya no son (reuse1, ffr D=3), que ' ...
     'son los que ejecuto T2. El ancla compararia filas distintas.']);
fprintf('T2 publicado: step_km = %g, Nt = %s\n', PUB.T2.step_km, mat2str(PUB.T2.Nt));
fprintf('NOTA: el original eligio nChunks = %s por RAM libre; aqui se fuerza %d.\n', ...
    mat2str(PUB.T2.nChunks), FORCE_CHUNK);
fprintf('      Si el troceado es exacto, el resultado debe ser IDENTICO igualmente.\n');
fprintf(['NOTA: el barrido corre con %d esquemas y el ancla compara SOLO las filas\n' ...
         '      1:%d (las publicadas). Que salga 0 prueba ademas que anadir\n' ...
         '      ffr(D=4) no perturba a los dos esquemas anteriores.\n\n'], nSCH, nPUB);

A6 = run_nt_sweep(cfg0, SCH, ST, STEP_ANCLA, FORCE_CHUNK);

print_p5_table(sprintf('ANCLA (step_km = %d)', STEP_ANCLA), ST, SCH, A6, false);

dAnc = max(max(abs(A6.p5(1:nPUB,:) - PUB.T2.p5)));
fprintf('\nmax|dif| vs T2 publicado (filas 1:%d) = %.3e dB\n', nPUB, dAnc);
if dAnc ~= 0
    for s = 1:nPUB
        fprintf('\n  %-14s obtenido = %s\n  %-14s publicado= %s\n', ...
            SCH{s}.name, mat2str(round(A6.p5(s,:),6)), '', ...
            mat2str(round(PUB.T2.p5(s,:),6)));
    end
    error('run_convergence_nt_step4:ancla', ...
        ['El ancla NO reproduce T2 (max|dif| = %.6g dB). O la copia de eval_case ' ...
         'ha divergido de la original, o el troceado no es exacto. PARAR: la ' ...
         'comparacion de la PARTE 2 no seria interpretable.'], dAnc);
end
fprintf('ANCLA OK: la copia de eval_case reproduce T2 EXACTAMENTE (max|dif| = 0).\n');

%% ------------------------------------------------------------------------
%  PARTE 2 -- Barrido de Nt con la rejilla de TRABAJO
%  ------------------------------------------------------------------------
fprintf('\n=========== PARTE 2: barrido de Nt con step_km = %d ===========\n', STEP_TRABAJO);
A4 = run_nt_sweep(cfg0, SCH, ST, STEP_TRABAJO, FORCE_CHUNK);

print_p5_table(sprintf('REJILLA DE TRABAJO (step_km = %d)', STEP_TRABAJO), ...
    ST, SCH, A4, true);
fprintf('(pico[GB] = pico modelado de una pasada ENTERA sin trocear; el troceado lo acota)\n');

% Tabla de GEOMETRIA VALIDA (elev >= 45 deg), informativa, en su propio bloque para
% no ensanchar la principal ahora que hay tres esquemas.
fprintf('\n-- misma tabla restringida a geometria valida (elev >= %d deg) --\n', ST.elMin_valid);
fprintf('%-6s', 'Nt');
for s = 1:nSCH, fprintf(' %16s', SCH{s}.name); end
fprintf('\n%s\n', repmat('-', 1, 6 + 17*nSCH));
for i = 1:numel(ST.Nt_list)
    fprintf('%-6d', ST.Nt_list(i));
    for s = 1:nSCH, fprintf(' %+16.4f', A4.p5V(s,i)); end
    fprintf('\n');
end

%% ------------------------------------------------------------------------
%  PARTE 3 -- Criterio de convergencia (EL YA EXISTENTE, sin tocar)
%  ------------------------------------------------------------------------
fprintf('\n=========== PARTE 3: criterio de convergencia (tol = %.1f dB) ===========\n', ST.tol_dB);
fprintf('Regla (converged_step, copia literal): el Nt de MENOR coste cuyo p5 difiere\n');
fprintf('menos de %.1f dB de TODOS los Nt mayores, tomando el maximo de los esquemas.\n', ST.tol_dB);
fprintf('Se reporta el agregado sobre los %d esquemas PUBLICADOS (cifra historica,\n', nPUB);
fprintf('trazable) y sobre los %d actuales, que es el que decide ahora.\n\n', nSCH);

R = struct();
for w = 1:2
    if w == 1, A = A6; lab = sprintf('step_km = %d (ANCLA, = T2 publicado)', STEP_ANCLA);
    else,      A = A4; lab = sprintf('step_km = %d (REJILLA DE TRABAJO)', STEP_TRABAJO); end

    fprintf('--- %s ---\n', lab);
    fprintf('%-6s', 'Nt');
    for s = 1:nSCH, fprintf(' %14s', ['max|d| ' SCH{s}.short]); end
    fprintf(' %10s %10s %8s\n', sprintf('max %d pub', nPUB), sprintf('max %d', nSCH), 'cumple');
    fprintf('%s\n', repmat('-', 1, 6 + 15*nSCH + 32));
    for i = 1:numel(ST.Nt_list)
        fprintf('%-6d', ST.Nt_list(i));
        if i < numel(ST.Nt_list)
            d = nan(1,nSCH);
            for s = 1:nSCH
                d(s) = max(abs(A.p5(s,i+1:end) - A.p5(s,i)));
                fprintf(' %14.4f', d(s));
            end
            dpub = max(d(1:nPUB));  dall = max(d);
            fprintf(' %10.4f %10.4f %8s\n', dpub, dall, yesno(dall < ST.tol_dB));
        else
            for s = 1:nSCH, fprintf(' %14s', '-'); end
            fprintf(' %10s %10s %8s\n', '-','-','-');
        end
    end

    [ntcP, dp5P] = converged_step(ST.Nt_list, A.p5(1:nPUB,:), ST.tol_dB, A.ok);
    [ntc,  dp5 ] = converged_step(ST.Nt_list, A.p5,           ST.tol_dB, A.ok);
    fprintf('difs consecutivas (max %d publicados): %s\n', nPUB, mat2str(round(dp5P,4)));
    fprintf('difs consecutivas (max %d esquemas):   %s\n', nSCH, mat2str(round(dp5,4)));
    fprintf('>> Nt CONVERGIDO con %d esquemas: %s\n', nPUB, conv_txt(ntcP, ST));
    fprintf('>> Nt CONVERGIDO con %d esquemas: %s\n\n', nSCH, conv_txt(ntc, ST));

    if w == 1
        R.Nt_conv_step6 = ntc;  R.dp5_step6 = dp5;
        R.Nt_conv_step6_pub = ntcP;  R.dp5_step6_pub = dp5P;
    else
        R.Nt_conv_step4 = ntc;  R.dp5_step4 = dp5;
        R.Nt_conv_step4_pub = ntcP;  R.dp5_step4_pub = dp5P;
    end
end

% La pregunta CONCRETA: ¿sigue Nt = 61 cumpliendo el criterio a step_km = 4?
i61 = find(ST.Nt_list == 61, 1);
d61_step4 = max(max(abs(A4.p5(:,i61+1:end) - A4.p5(:,i61))));
d61_step6 = max(max(abs(A6.p5(:,i61+1:end) - A6.p5(:,i61))));
R.d61_step4 = d61_step4;  R.d61_step6 = d61_step6;
R.d61_step4_sch = max(abs(A4.p5(:,i61+1:end) - A4.p5(:,i61)), [], 2)';   % por esquema
R.d61_step6_sch = max(abs(A6.p5(:,i61+1:end) - A6.p5(:,i61)), [], 2)';
fprintf('PREGUNTA DIRECTA -- ¿Nt = 61 cumple el criterio con la rejilla de trabajo?\n');
fprintf('  max|p5(Nt) - p5(61)| para Nt > 61, POR ESQUEMA (step_km = %d):\n', STEP_TRABAJO);
for s = 1:nSCH
    fprintf('    %-16s %.4f dB  ->  %s\n', SCH{s}.name, R.d61_step4_sch(s), ...
        yesno_txt(R.d61_step4_sch(s) < ST.tol_dB, 'CUMPLE', 'NO CUMPLE'));
end
fprintf('  agregado:  step_km=%d -> %.4f dB | step_km=%d -> %.4f dB  (umbral %.1f dB)\n', ...
    STEP_ANCLA, d61_step6, STEP_TRABAJO, d61_step4, ST.tol_dB);
fprintf('  VEREDICTO a step_km=%d sobre los %d esquemas: %s\n', STEP_TRABAJO, nSCH, ...
    yesno_txt(d61_step4 < ST.tol_dB, 'CUMPLE', 'NO CUMPLE'));

%% ------------------------------------------------------------------------
%  PARTE 4 -- Cuanto se mueve cada punto al cambiar de rejilla
%  ------------------------------------------------------------------------
fprintf('\n=========== PARTE 4: step_km %d -> %d, punto a punto ===========\n', ...
    STEP_ANCLA, STEP_TRABAJO);
dP5 = A4.p5 - A6.p5;              % [nSCH x nNt] desplazamiento por cambio de rejilla
dR1 = dP5(1,:);  dFF = dP5(2,:);  % nombres historicos (se siguen guardando)

fprintf('%-6s', 'Nt');
for s = 1:nSCH, fprintf(' %11s %11s %10s', [SCH{s}.short ' step6'], ...
        [SCH{s}.short ' step4'], ['d ' SCH{s}.short]); end
fprintf('\n%s\n', repmat('-', 1, 6 + 35*nSCH));
for i = 1:numel(ST.Nt_list)
    fprintf('%-6d', ST.Nt_list(i));
    for s = 1:nSCH
        fprintf(' %11.4f %11.4f %+10.4f', A6.p5(s,i), A4.p5(s,i), dP5(s,i));
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('-', 1, 6 + 35*nSCH));
for s = 1:nSCH
    [~, imax] = max(abs(dP5(s,:)));
    fprintf('%-16s desplazamiento medio %+.4f dB | maximo %+.4f dB | recorrido del eje Nt %.4f dB\n', ...
        SCH{s}.name, mean(dP5(s,:)), dP5(s,imax), ...
        max(A4.p5(s,:)) - min(A4.p5(s,:)));
end
fprintf(['(el "recorrido del eje Nt" es a step_km=%d; comparar con el desplazamiento\n' ...
         ' de rejilla para ver si encadenar los ejes es gratis)\n'], STEP_TRABAJO);

%% ------------------------------------------------------------------------
%  PARTE 5 -- Los valores PUBLICADOS a step_km = 4 no se pueden haber movido
%  ------------------------------------------------------------------------
%  Segunda puerta, independiente del ancla: el ancla comprueba la rejilla de 6 km
%  contra el .mat; esto comprueba la de 4 km contra las cifras publicadas del
%  estudio. Si anadir ffr(D=4) hubiera perturbado algo, cae aqui.
fprintf('\n=========== PARTE 5: los %d esquemas publicados no se han movido ===========\n', nPUB);
fprintf('%-6s', 'Nt');
for s = 1:nPUB, fprintf(' %11s %11s %10s', [SCH{s}.short ' pub'], ...
        [SCH{s}.short ' hoy'], ['dif ' SCH{s}.short]); end
fprintf('\n%s\n', repmat('-', 1, 6 + 35*nPUB));
for i = 1:numel(ST.Nt_list)
    fprintf('%-6d', ST.Nt_list(i));
    for s = 1:nPUB
        fprintf(' %11.4f %11.4f %10.2e', PUB_STEP4(s,i), A4.p5(s,i), ...
            abs(A4.p5(s,i) - PUB_STEP4(s,i)));
    end
    fprintf('\n');
end
dPub = max(max(abs(A4.p5(1:nPUB,:) - PUB_STEP4)));
fprintf('%s\n', repmat('-', 1, 6 + 35*nPUB));
fprintf('max|dif| vs CLAUDE.md = %.3e dB (tolerancia %.0e = redondeo a 4 decimales)\n', ...
    dPub, PUB_TOL);
if ~(dPub <= PUB_TOL)
    error('run_convergence_nt_step4:publicados', ...
        ['Los valores publicados de reuse1 / ffr(D=3) a step_km = %d SE HAN MOVIDO ' ...
         '(max|dif| = %.6g dB > %.0e). PARAR: anadir el tercer esquema no puede ' ...
         'cambiarlos.'], STEP_TRABAJO, dPub, PUB_TOL);
end
fprintf('OK: reuse1 y ffr(D=3) reproducen lo publicado dentro del redondeo.\n');

%% ------------------------------------------------------------------------
%  PARTE 6 -- Figura del barrido (las tres curvas)
%  ------------------------------------------------------------------------
%  TRES paneles porque la escala manda: los esquemas estan separados ~12 dB y en
%  valor absoluto la convergencia (decimas de dB) no se ve. Los paneles (a) y (b)
%  dan el valor absoluto en las dos rejillas; el (c) da la MAGNITUD QUE DECIDE, que
%  es la del criterio: max|p5(Nt) - p5(Nt')| sobre TODOS los Nt' mayores.
FIGDIR = 'figs_faseA';
mk = {'o-','s-','^-'};
fF = figure('Name','Convergencia vs Nt (3 esquemas)','Color','w','Visible','off', ...
            'Position',[100 100 1420 420]);

for w = 1:2
    if w == 1, A = A4; stp = STEP_TRABAJO; pan = '(a) rejilla de TRABAJO';
    else,      A = A6; stp = STEP_ANCLA;   pan = '(b) ANCLA';  end
    ax = subplot(1,3,w);  hold(ax,'on');
    for s = 1:nSCH
        plot(ax, ST.Nt_list(A.ok), A.p5(s,A.ok), mk{min(s,end)}, ...
            'LineWidth',1.6, 'DisplayName', SCH{s}.name);
    end
    xline(ax, 61, 'k:', 'LineWidth',1.4, 'HandleVisibility','off');
    grid(ax,'on'); box(ax,'on');
    set(ax,'XTick',ST.Nt_list); xlim(ax,[ST.Nt_list(1)-12, ST.Nt_list(end)+12]);
    xlabel(ax,'numero de instantes N_t'); ylabel(ax,'SINR_{edge,p5} [dB]');
    title(ax, sprintf('%s: step\\_km = %d (M = %d)', pan, stp, A.M(1)));
    legend(ax,'Location','east');
end

% (c) la magnitud del criterio, a la rejilla de trabajo
dcrit = nan(nSCH, numel(ST.Nt_list));
for i = 1:numel(ST.Nt_list)-1
    dcrit(:,i) = max(abs(A4.p5(:,i+1:end) - A4.p5(:,i)), [], 2);
end
ax3 = subplot(1,3,3);  hold(ax3,'on');
for s = 1:nSCH
    plot(ax3, ST.Nt_list(1:end-1), dcrit(s,1:end-1), mk{min(s,end)}, ...
        'LineWidth',1.6, 'DisplayName', SCH{s}.name);
end
yline(ax3, ST.tol_dB, 'k--', 'LineWidth',1.6, ...
      'Label',sprintf('criterio = %.1f dB', ST.tol_dB), 'FontSize',10, ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
xline(ax3, 61, 'k:', 'LineWidth',1.4, 'HandleVisibility','off');
grid(ax3,'on'); box(ax3,'on');
set(ax3,'XTick',ST.Nt_list(1:end-1));
xlim(ax3,[ST.Nt_list(1)-12, ST.Nt_list(end-1)+12]);
ylim(ax3,[0, max(0.34, 1.15*max(dcrit(:)))]);
xlabel(ax3,'numero de instantes N_t'); ylabel(ax3,'max_{N_t'' > N_t} |\Delta p5| [dB]');
title(ax3, sprintf('(c) criterio de convergencia (step\\_km = %d)', STEP_TRABAJO));
legend(ax3,'Location','northeast');

fpng = save_fig(fF, FIGDIR, 'faseA_d_convergencia_vs_Nt_step4_3esquemas');
close(fF);
fprintf('\n[figura] %s\n', fpng);

%% 5. Guardar
timing = struct('warmup_p618_s', tWarm, 'study_s', toc(tStudy));
fprintf('\n[tiempos] barrido %.1f s (%.1f min) (+ %.1f s de warmup P.618)\n', ...
    timing.study_s, timing.study_s/60, timing.warmup_p618_s);
fprintf('[memoria] pico del proceso: %.2f GB | RAM libre ahora: %.2f GB\n', ...
    proc_peak_GB(), ram_free_GB());

save('convergence_nt_step4.mat', 'ST','cfg0','SCH','A4','A6','R', ...
     'STEP_ANCLA','STEP_TRABAJO','FORCE_CHUNK','dAnc','dR1','dFF','dP5', ...
     'nSCH','nPUB','PUB_STEP4','dPub','timing');
fprintf('\nResultados guardados en convergence_nt_step4.mat\n');
fprintf('(convergence_studyA.mat NO se ha tocado)\n');

% =========================================================================
% ============================ FUNCIONES ==================================
% =========================================================================

function A = run_nt_sweep(cfg0, SCH, ST, step_km, forceChunk)
%RUN_NT_SWEEP  Barrido de Nt a un step_km dado. Solo orquesta: la evaluacion la
%   hace eval_case, copia literal de run_convergence_studyA.
n  = numel(ST.Nt_list);
nS = numel(SCH);                       % 3 desde la ampliacion; eval_case ya era generico
A = struct('Nt',ST.Nt_list, 'step_km',step_km, 'dt',nan(1,n), 'M',nan(1,n), ...
           'Medge',nan(1,n), 'p5',nan(nS,n), 'p5V',nan(nS,n), ...
           'tens_GB',nan(1,n), 'peak_model_GB',nan(1,n), 'nChunks',nan(1,n), ...
           'ok',false(1,n), 'secs',nan(1,n));
for i = 1:n
    cfg = cfg0;
    cfg.ground.step_km = step_km;
    cfg.time.dt        = ST.window_s / (ST.Nt_list(i) - 1);

    E = eval_case(cfg, SCH, ST, forceChunk);
    A.dt(i)            = cfg.time.dt;
    A.M(i)             = E.M;
    A.tens_GB(i)       = E.tens_GB;
    A.peak_model_GB(i) = E.peak_model_GB;
    A.nChunks(i)       = E.nChunks;
    A.secs(i)          = E.secs;
    A.ok(i)            = E.ok;
    if ~E.ok
        fprintf('  Nt=%-4d SALTADO: %s\n', ST.Nt_list(i), E.why);
        continue;
    end
    A.Medge(i) = E.K{1}.viab.n_edge;
    for s = 1:nS
        A.p5(s,i)  = E.K{s}.viab.SINR_edge_p5;
        A.p5V(s,i) = E.KV{s}.viab.SINR_edge_p5;
    end
    fprintf('  step=%d Nt=%-4d M=%-5d ->', step_km, ST.Nt_list(i), E.M);
    for s = 1:nS, fprintf(' %s %+8.4f |', SCH{s}.name, A.p5(s,i)); end
    fprintf('  (%d bloques, %.1f s)\n', E.nChunks, E.secs);
end
end

% -------------------------------------------------------------------------
function print_p5_table(lab, ST, SCH, A, withMem)
%PRINT_P5_TABLE  Tabla de p5 por Nt con UNA columna por esquema. Solo imprime.
%   Se escribe una vez y la usan las PARTES 1 y 2, para que las dos tablas no
%   puedan divergir al cambiar el numero de esquemas.
nS = numel(SCH);
fprintf('\n-- %s --\n', lab);
fprintf('%-6s %8s %8s %10s', 'Nt','dt[s]','M','M_borde');
for s = 1:nS, fprintf(' %16s', SCH{s}.name); end
if withMem, fprintf(' %8s %9s %8s', 'bloques','pico[GB]','t[s]'); end
fprintf('\n%s\n', repmat('-', 1, 34 + 17*nS + 28*double(withMem)));
for i = 1:numel(ST.Nt_list)
    fprintf('%-6d %8.0f %8d %10d', ST.Nt_list(i), A.dt(i), A.M(i), A.Medge(i));
    for s = 1:nS, fprintf(' %+16.4f', A.p5(s,i)); end
    if withMem
        fprintf(' %8d %9.2f %8.1f', A.nChunks(i), A.peak_model_GB(i), A.secs(i));
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('-', 1, 34 + 17*nS + 28*double(withMem)));
end

% -------------------------------------------------------------------------
function s = conv_txt(ntc, ST)
if isnan(ntc)
    s = sprintf('NO converge en el rango barrido con tol = %.1f dB', ST.tol_dB);
else
    s = sprintf('Nt = %d (dt = %.0f s)', ntc, ST.window_s/(ntc-1));
end
end

% -------------------------------------------------------------------------
% A PARTIR DE AQUI: COPIA LITERAL de las funciones locales de
% run_convergence_studyA.m. No modificar sin modificar alli tambien; el ancla
% de la PARTE 1 detiene el script si divergen.
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
    error('run_convergence_nt_step4:memModel', ...
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
