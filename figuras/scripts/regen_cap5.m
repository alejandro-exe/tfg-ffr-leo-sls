%% REGEN_CAP5   Figuras 5.1 a 5.6, SIN TITULO, a REDACCION/_staging y promocion.
%
%  Reglas aplicadas: sin title()/sgtitle(); paneles etiquetados (a)/(b)/(c); ejes con
%  unidad; texto pintado en ASCII; fuentes >= 11 pt; 300 dpi via save_fig.
%
%  REGLA DE RECALCULO: ninguno de estos bloques llama a compute_sinr_ffr,
%  compute_interference, compute_kpis, ffr_allocate, ffr_policy, run_one_density ni
%  run_sweep_points. Todo lo graficado sale del .mat, salvo lo declarado en §5.2 y
%  §5.4 (re-derivacion GEOMETRICA / patron analitico).
%
%  Uso: run('.../figuras/scripts/regen_cap5.m')

HERE     = fileparts(mfilename('fullpath'));
PROJROOT = fileparts(fileparts(HERE));
addpath(PROJROOT);
STAGE  = fullfile(PROJROOT,'REDACCION','_staging');
FINAL  = fullfile(PROJROOT,'REDACCION','figuras');
if ~exist(STAGE,'dir'), mkdir(STAGE); end

promote = @(name) copyfile(fullfile(STAGE,[name '.png']), fullfile(FINAL,[name '.png']));
bigfont = @(fh) set(findall(fh,'-property','FontSize'),'FontSize',11);
panel   = @(ax,txt) text(ax, 0.015, 0.985, txt, 'Units','normalized', ...
             'FontSize',12,'FontWeight','bold','VerticalAlignment','top', ...
             'HorizontalAlignment','left','BackgroundColor',[1 1 1 ],'Margin',1);

fprintf('\n################ CAPITULO 5 ################\n');

%% ===================== 5.1  Calibracion vs 3GPP =========================
%  Origen: export_e0_figures.m:111-128 (fig7_calib_vs_3gpp). Categoria A.
E0 = load(fullfile(PROJROOT,'calibration_e0_results.mat'), 'CIR');

% Constantes de referencia 3GPP TR 38.821 V16.0.0, Tabla 6.1.1.2-1, transcritas
% LITERALMENTE de export_e0_figures.m:36 (estan hardcodeadas alli, no en el .mat).
ref_SIR = [-3.0, -1.0,  1.2];                       % [p5 p50 p95]
own_SIR = [E0.CIR.stats.CIR_p5, E0.CIR.stats.CIR_p50, E0.CIR.stats.CIR_p95];
dSIR    = own_SIR - ref_SIR;

fprintf('\n--- 5.1 ---\n');
fprintf('  ref 3GPP SIR  [p5 p50 p95] = [%.4f %.4f %.4f]\n', ref_SIR);
fprintf('  propio CIR    [p5 p50 p95] = [%.4f %.4f %.4f]\n', own_SIR);
fprintf('  delta                      = [%+.4f %+.4f %+.4f]  max|d| = %.4f dB\n', ...
        dSIR, max(abs(dSIR)));

% Verificacion NUMERICA contra los valores PUBLICADOS del corpus
% ("CIR media -0.53, p5 -2.28, p50 -0.59, p95 +1.28"; "max|D| vs 3GPP = 0.72 dB")
ref_pub = [-2.28, -0.59, 1.28];
d_pub   = max(abs(own_SIR - ref_pub));
fprintf('  vs valores publicados [-2.28 -0.59 +1.28]: max|dif| = %.4f dB\n', d_pub);
assert(d_pub <= 0.005, '5.1: los CIR del .mat no reproducen los publicados (%.4f dB).', d_pub);
d_max = max(abs(dSIR));
assert(abs(d_max - 0.72) <= 0.005, '5.1: max|delta| vs 3GPP = %.4f, publicado 0.72.', d_max);
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off','Position',[100 100 900 420]);
ax1 = subplot(1,2,1);
b = bar([ref_SIR; own_SIR]', 'grouped');
b(1).FaceColor=[0.85 0.7 0.2]; b(2).FaceColor=[0.2 0.4 0.85];
set(gca,'XTickLabel',{'p5','p50','p95'}); grid on; box on;
xlabel('percentil'); ylabel('DL Geometry SIR / CIR [dB]');
% leyenda al SURESTE: al noroeste tapaba la etiqueta de panel (a)
legend({'3GPP TR 38.821','propio (CIR)'},'Location','southeast');
panel(ax1,'(a)');
ax2 = subplot(1,2,2);
bd = bar(dSIR); bd.FaceColor=[0.3 0.6 0.3];
set(gca,'XTickLabel',{'p5','p50','p95'}); grid on; box on;
xlabel('percentil'); ylabel('propio - 3GPP [dB]'); ylim([0 1.2]);
yline(1,'r--','tol ~1 dB');
for i=1:3
    text(i, dSIR(i)+0.05, sprintf('%+.2f',dSIR(i)),'HorizontalAlignment','center','FontSize',11);
end
panel(ax2,'(b)');
bigfont(f); save_fig(f, STAGE, 'fig_5_01_calibracion_vs_3gpp'); close(f);
promote('fig_5_01_calibracion_vs_3gpp');

%% ===================== 5.2  Validacion geometrica ========================
%  Origen: export_e0_figures.m:130-148 (fig6_geom_validation). Categoria B:
%  el .mat NO guarda V; se RE-DERIVA con validate_geometry_satscenario(cfg), que
%  produce ELEVACION y RANGO -- magnitudes puramente GEOMETRICAS, sin ningun KPI.
%  La fidelidad se verifica contra los residuos PUBLICADOS del corpus.
E0c = load(fullfile(PROJROOT,'calibration_e0_results.mat'), 'cfg');
cfgE0 = E0c.cfg;
tvecE0 = cfgE0.time.t0 : cfgE0.time.dt : cfgE0.time.t0 + cfgE0.time.duration;

fprintf('\n--- 5.2 ---\n');
okV = true;
try
    V = validate_geometry_satscenario(cfgE0);
    close all;                       % la validacion abre su propia figura; se descarta
catch ME
    okV = false;
    fprintf('  ** NO se pudo re-derivar V: %s\n', ME.message);
end

if okV
    % VALORES PUBLICADOS: residuos sobre jsat=322, el servidor en su maxima
    % elevacion (paso REAL: elevacion [-87.37, +88.83] deg, muy por encima de la
    % mascara de 30 deg).
    %
    % Los valores anteriores (max|dEl| 0.037 / rms 0.023 deg; max|dRng| 6.0 /
    % rms 3.4 km) correspondian a jsat=1064, que NUNCA estaba sobre el horizonte
    % (elevacion [-58.92, -19.42] deg en toda la ventana): eran un ARTEFACTO del
    % bug de seleccion ya corregido en validate_geometry_satscenario (jsat se
    % elegia con theta_g0 por defecto y luego se trazaba con la epoca alineada
    % theta_g0 = GMST0, dos convenciones distintas). No hubo cambio de modelo ni
    % de metodo de comparacion.
    fprintf('  max|dEl|  = %.4f deg (publicado 0.2815) | rms = %.4f deg (publicado 0.0178)\n', ...
            V.maxAbsEl_deg, V.rmsEl_deg);
    fprintf('  max|dRng| = %.3f km  (publicado 1.284)  | rms = %.3f km  (publicado 0.623)\n', ...
            V.maxAbsRng_km, V.rmsRng_km);
    dEl_pub  = abs(V.maxAbsEl_deg - 0.2815);
    dRng_pub = abs(V.maxAbsRng_km - 1.284);
    fprintf('  desviacion vs publicado: elevacion %.4f deg, rango %.3f km\n', dEl_pub, dRng_pub);
    if dEl_pub > 0.001 || dRng_pub > 0.05
        okV = false;
        fprintf('  ** La re-derivacion NO reproduce los residuos publicados. NO se promueve.\n');
    else
        fprintf('  VERIFICACION OK\n');
    end
end

if okV
    tmin = tvecE0/60;
    f = figure('Color','w','Visible','off','Position',[100 100 760 720]);
    ax1 = subplot(3,1,1);
    plot(tmin,V.el_own,'b-','LineWidth',1.5); hold on;
    plot(tmin,V.el_tb,'r--','LineWidth',1.2);
    grid on; box on; ylabel('elevacion [deg]');
    legend({'propia','satelliteScenario'},'Location','best'); panel(ax1,'(a)');
    ax2 = subplot(3,1,2);
    plot(tmin,V.rng_own,'b-','LineWidth',1.5); hold on;
    plot(tmin,V.rng_tb,'r--','LineWidth',1.2);
    grid on; box on; ylabel('rango [km]');
    legend({'propia','satelliteScenario'},'Location','best'); panel(ax2,'(b)');
    ax3 = subplot(3,1,3);
    yyaxis left
    plot(tmin, V.el_own-V.el_tb, '-', 'LineWidth',1.5);
    ylabel('\Delta elevacion [deg]');
    yyaxis right
    plot(tmin, V.rng_own-V.rng_tb, '--', 'LineWidth',1.5);
    ylabel('\Delta rango [km]');
    grid on; box on; xlabel('tiempo [min]');
    legend({'\Delta elevacion [deg]','\Delta rango [km]'},'Location','best');
    panel(ax3,'(c)');
    bigfont(f); save_fig(f, STAGE, 'fig_5_02_validacion_geometria'); close(f);
    promote('fig_5_02_validacion_geometria');
else
    fprintf('  -> se CONSERVA la copia con titulo ya presente en figuras/.\n');
end

%% ===================== 5.3  Convergencia vs rejilla ======================
%  Origen: run_convergence_studyA.m:552-559. Categoria A.
FA = load(fullfile(PROJROOT,'convergence_studyA.mat'), 'ST','T1','T4');

fprintf('\n--- 5.3 ---\n');
fprintf('  steps_km = %s\n', mat2str(FA.ST.steps_km));
fprintf('  reuse1   = %s\n', mat2str(round(FA.T1.p5(1,FA.T1.ok),3)));
fprintf('  ffr(D=3) = %s\n', mat2str(round(FA.T1.p5(2,FA.T1.ok),3)));
% Verificacion contra la tabla publicada del estudio de convergencia:
%   step 8 -> -4.997/+3.645 · 6 -> -5.038/+3.757 · 4 -> -4.950/+3.865
%   3 -> -4.972/+3.782 · 2 -> -4.974/+3.812
pub53 = [-4.997 -5.038 -4.950 -4.972 -4.974;   % reuse1, en orden de steps_km
          3.645  3.757  3.865  3.782  3.812];  % ffr
[sk, iso] = sort(FA.ST.steps_km, 'descend');
got53 = [FA.T1.p5(1,iso); FA.T1.p5(2,iso)];
d53 = max(abs(got53(:) - pub53(:)));
fprintf('  steps ordenados %s -> max|dif| vs publicado = %.4f dB\n', mat2str(sk), d53);
assert(d53 <= 0.0006, '5.3: no reproduce la tabla de FASE A (max|dif| = %.4f).', d53);
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off');
plot(FA.ST.steps_km(FA.T1.ok), FA.T1.p5(1,FA.T1.ok), 'o-','LineWidth',1.6); hold on;
plot(FA.ST.steps_km(FA.T1.ok), FA.T1.p5(2,FA.T1.ok), 's-','LineWidth',1.6);
set(gca,'XDir','reverse'); grid on; box on;
xlabel('paso de la rejilla [km] (mas FINA hacia la derecha)');
ylabel('SINR_{edge,p5} [dB]');
legend({'reuse1','ffr (\Delta=3)'},'Location','best');
bigfont(f); save_fig(f, STAGE, 'fig_5_03_convergencia_rejilla'); close(f);
promote('fig_5_03_convergencia_rejilla');

%% ===================== 5.4  Suelo de lobulos vs Bessel ===================
%  Origen: run_convergence_nrings.m:642-663. Categoria A: los datos del .mat son
%  cfg y RES; las DOS CURVAS las evalua beam_gain_dB en tiempo de graficado, igual
%  que hace el original (patron analitico, no es un indicador).
NR = load(fullfile(PROJROOT,'convergence_nrings.mat'), 'RES','cfg');
sdeg  = NR.RES{end}.spacing_deg;
HPBWd = NR.cfg.radio.beamwidth3dB_deg;
flr   = NR.cfg.radio.sidelobe_floor_dB;
rr    = 0:NR.RES{end}.nRings;
gBes  = arrayfun(@(r) beam_gain_dB(r*sdeg, HPBWd, -Inf), rr);
gUsed = arrayfun(@(r) beam_gain_dB(r*sdeg, HPBWd, flr),  rr);
iOn   = find(gBes < flr, 1);

fprintf('\n--- 5.4 ---\n');
fprintf('  s = %.4f deg | HPBW = %.4f deg | suelo = %g dB | anillos 0..%d\n', ...
        sdeg, HPBWd, flr, rr(end));
fprintf('  Bessel puro por anillo = %s\n', mat2str(round(gBes,2)));
fprintf('  el suelo MANDA desde el anillo %d\n', rr(iOn));
% Publicado: "El Bessel real vale -28.91 dB en el anillo 4 pero -35.03 / -45.56 /
% -49.71 dB en los anillos 5 / 6 / 7 -> el suelo entra exactamente en el anillo 5"
pub54 = [-28.91 -35.03 -45.56 -49.71];      % anillos 4,5,6,7
d54 = max(abs(gBes(5:8) - pub54));
fprintf('  anillos 4..7 vs publicado: max|dif| = %.4f dB\n', d54);
assert(d54 <= 0.006, '5.4: el patron no reproduce los valores publicados (%.4f dB).', d54);
assert(rr(iOn) == 5, '5.4: el suelo manda desde el anillo %d, publicado 5.', rr(iOn));
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off');
plot(rr, gBes,  'o-', 'LineWidth',1.8); hold on;
plot(rr, gUsed, 's--','LineWidth',1.8);
yline(flr, 'k:', 'LineWidth',1.6, 'Label',sprintf('suelo de lobulos = %g dB', flr), ...
      'FontSize',11, 'LabelHorizontalAlignment','left');
xline(rr(iOn), 'r--', 'LineWidth',1.8, ...
      'Label',sprintf('el suelo MANDA desde el anillo %d', rr(iOn)), ...
      'FontSize',11, 'LabelVerticalAlignment','bottom');
grid on; box on;
xlabel('anillo co-canal r'); ylabel('ganancia relativa del haz co-canal [dB]');
legend({'Bessel real (sin suelo)','patron USADO (con suelo)'},'Location','southwest');
bigfont(f); save_fig(f, STAGE, 'fig_5_04_suelo_lobulos_vs_bessel'); close(f);
promote('fig_5_04_suelo_lobulos_vs_bessel');

%% ============ 5.5  Convergencia vs numero de instantes ===================
%  Origen: run_convergence_studyA.m:561-567 (figura faseA_b_convergencia_vs_
%  numero_instantes.png), PROMOVIDA del anexo B al cuerpo del capitulo.
%  Categoria A: todo sale del .mat, sin recalcular nada.
%
%  AVISO DE PROCEDENCIA (hay que decirlo en el pie de figura): la TAREA 2 de la
%  FASE A se ejecuto con el step CONVERGIDO de la tarea 1, que es
%  step_km = 6 (M = 137), no con el valor de trabajo de 4 km
%  (run_convergence_studyA.m:215, step_use = step_conv, y el .mat guarda
%  T2.step_km = 6). Se comprueba abajo para que no se cite mal.
FB = load(fullfile(PROJROOT,'convergence_studyA.mat'), 'T2','ST');

fprintf('\n--- 5.5 ---\n');
fprintf('  step_km del barrido de Nt = %g (NO es el valor de trabajo de 4 km)\n', FB.T2.step_km);
fprintf('  Nt       = %s\n', mat2str(FB.T2.Nt));
fprintf('  reuse1   = %s\n', mat2str(round(FB.T2.p5(1,:),4)));
fprintf('  ffr(D=3) = %s\n', mat2str(round(FB.T2.p5(2,:),4)));

% Constantes de referencia: tabla de instantes de FASE A (convergence_studyA.mat,
% T2), ya verificada contra el .mat. NO se ajustan para que pasen los assert.
pub55_Nt = [31 61 121 241];
pub55    = [-5.0928 -5.0378 -5.0317 -4.9973;    % reuse1, en orden de Nt creciente
             3.6073  3.7570  3.7510  3.8538];   % ffr (Delta = 3, alpha = 0.4)

assert(all(FB.T2.ok), '5.5: hay puntos NO validos en T2.ok = %s.', mat2str(FB.T2.ok));
[nt55, iso55] = sort(FB.T2.Nt, 'ascend');
assert(isequal(nt55, pub55_Nt), ...
    '5.5: los Nt del .mat (%s) no son los de referencia (%s).', ...
    mat2str(nt55), mat2str(pub55_Nt));
got55 = [FB.T2.p5(1,iso55); FB.T2.p5(2,iso55)];
d55   = max(abs(got55(:) - pub55(:)));
fprintf('  Nt ordenados %s -> max|dif| vs referencia = %.4f dB\n', mat2str(nt55), d55);
assert(d55 <= 0.001, ...
    ['5.5: no reproduce la tabla de instantes de FASE A (max|dif| = %.4f dB > 0.001).\n' ...
     '     obtenido reuse1 = %s\n     esperado  reuse1 = %s\n' ...
     '     obtenido ffr    = %s\n     esperado  ffr    = %s'], ...
    d55, mat2str(round(got55(1,:),4)), mat2str(pub55(1,:)), ...
    mat2str(round(got55(2,:),4)), mat2str(pub55(2,:)));
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off');
plot(nt55, got55(1,:), 'o-','LineWidth',1.6); hold on;
plot(nt55, got55(2,:), 's-','LineWidth',1.6);
grid on; box on;
% Eje X NO invertido: aqui "mas fino" es de forma natural mas a la derecha, al
% contrario que en la 5.3 (donde el eje es el PASO de la rejilla).
set(gca,'XTick', nt55);
xlim([nt55(1)-12, nt55(end)+12]);
xlabel('numero de instantes N_t');
ylabel('SINR_{edge,p5} [dB]');
legend({'reuse1','ffr (\Delta=3)'},'Location','best');
bigfont(f); save_fig(f, STAGE, 'fig_5_05_convergencia_instantes'); close(f);
promote('fig_5_05_convergencia_instantes');

%% ============ 5.6  Perfil de interferencia por anillo ====================
%  Origen: run_convergence_nrings.m:625-636 (figura nrings_b_perfil_
%  interferencia_por_anillo.png), PROMOVIDA del anexo B al cuerpo.
%  Categoria A: pct_ring sale del .mat; el acumulado es un cumsum de esa misma
%  serie, no un calculo nuevo.
%
%  PUNTO GRAFICADO: nRings = 7 (el cluster MAYOR) y esquema reuse1, que es lo que
%  pinta el original (RES{end}.scheme{1}). Aqui se localiza por VALOR y no por
%  'end'/indice fijo, para que la figura no cambie en silencio si algun dia se
%  reordena o se amplia el barrido.
%
%  EJE UNICO (no yyaxis) A PROPOSITO: el original no usa yyaxis y las dos series
%  estan en las MISMAS unidades (% de la potencia interferente total). Ponerlas en
%  dos ejes independientes permitiria leer mal la altura de las barras, asi que no
%  procede la etiqueta Y derecha. Las dos series se distinguen por ESTILO ademas de
%  por color (barra rellena vs linea DISCONTINUA con marcador hueco), mismo criterio
%  que la 5.2(c), para que la figura siga siendo legible impresa en blanco y negro.
NR6 = load(fullfile(PROJROOT,'convergence_nrings.mat'), 'RES');

iR7 = find(cellfun(@(R) R.nRings == 7, NR6.RES), 1);
assert(~isempty(iR7), '5.6: no hay ningun punto con nRings = 7 en convergence_nrings.mat.');
R7  = NR6.RES{iR7};
ic6 = find(cellfun(@(s) strcmp(s.name,'reuse1'), R7.scheme), 1);
assert(~isempty(ic6), '5.6: no se encuentra el esquema reuse1 en RES{%d}.scheme.', iR7);

pr6  = R7.scheme{ic6}.pct_ring(:).';
cum6 = cumsum(pr6);
rr6  = 0:numel(pr6)-1;

fprintf('\n--- 5.6 ---\n');
fprintf('  punto nRings = %d (%d haces), esquema %s\n', R7.nRings, R7.nBeams, R7.scheme{ic6}.name);
fprintf('  pct_ring  = %s\n', mat2str(round(pr6,3)));
fprintf('  acumulado = %s\n', mat2str(round(cum6,3)));

% Constantes de referencia: perfil por anillo de convergence_nrings.mat
% (RES{nRings=7}.scheme{reuse1}.pct_ring), ya verificado contra el .mat.
pub56     = [4.253 24.754 35.078 22.066 5.046 3.265 2.756 2.783];
pub56_cum = 100.000;

assert(numel(pr6) == numel(pub56), ...
    '5.6: el perfil tiene %d anillos y la referencia %d.', numel(pr6), numel(pub56));
d56  = max(abs(pr6 - pub56));
d56c = abs(cum6(end) - pub56_cum);
fprintf('  max|dif| vs referencia = %.4f %% | acumulado final %.4f %% (dif %.4f)\n', ...
        d56, cum6(end), d56c);
assert(d56 <= 0.01, ...
    ['5.6: el perfil por anillo no reproduce la referencia (max|dif| = %.4f %% > 0.01).\n' ...
     '     obtenido = %s\n     esperado = %s'], ...
    d56, mat2str(round(pr6,3)), mat2str(pub56));
assert(d56c <= 0.01, ...
    '5.6: el acumulado final es %.4f %% y la referencia %.3f %% (dif %.4f).', ...
    cum6(end), pub56_cum, d56c);
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off');
bar(rr6, pr6, 'FaceColor',[0.30 0.55 0.75], 'EdgeColor',[0.15 0.28 0.42]); hold on;
plot(rr6, cum6, 'ks--', 'LineWidth',1.8, 'MarkerSize',7, 'MarkerFaceColor','w');
yline(100, 'k:', 'LineWidth',1.2);
grid on; box on;
xlim([-0.6, rr6(end)+0.6]);
ylim([0 112]);
set(gca,'XTick', rr6);
xlabel('anillo co-canal r');
ylabel('aportacion a la interferencia [%]');
legend({'aportacion del anillo','acumulado'},'Location','east');
bigfont(f); save_fig(f, STAGE, 'fig_5_06_interferencia_por_anillo'); close(f);
promote('fig_5_06_interferencia_por_anillo');

fprintf('\n[cap 5] terminado.\n');
