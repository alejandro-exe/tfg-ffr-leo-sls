%% MAKE_FIG_4_03A_PATRON_HAZ   Candidata A a la figura 4.3
%  Corte del patron de haz del SATELITE (un solo haz), con dos curvas superpuestas:
%    - patron IDEAL : beam_gain_dB(theta, HPBW)                  (Bessel 38.811 puro)
%    - patron USADO : beam_gain_dB(theta, HPBW, sidelobe_floor)  (con suelo de lobulos)
%
%  PERFIL: config_default() (perfil Ku de ESTUDIO). NUNCA config_calib.
%
%  Que ilustra: el suelo plano de lobulos laterales es una DECISION de modelado
%  (cfg.radio.sidelobe_floor_dB = -30 dB) que acota la cola del Bessel ideal, cuyos
%  nulos profundos infraestiman la interferencia. Es el artefacto que documentan
%  CONVERGENCIA DE nRings (el suelo entra en el anillo 5) y APUNTAMIENTO DE LOS
%  INTERFERENTES (la rama 'nadir' opera enteramente sobre el suelo).
%
%  SIN TITULO (el pie va en LaTeX). Solo lee el motor; escribe unicamente en
%  REDACCION/figuras/.
%
%  Uso: desde la raiz del proyecto,
%       run('figuras/scripts/make_fig_4_03a_patron_haz.m')

% Rutas ancladas a la ubicacion del script (ver nota en regen_fig_4_02_*).
HERE     = fileparts(mfilename('fullpath'));
PROJROOT = fileparts(fileparts(HERE));
addpath(PROJROOT);
OUTDIR   = fullfile(PROJROOT,'REDACCION','figuras');

cfg   = config_default();
HPBW  = cfg.radio.beamwidth3dB_deg;
FLOOR = cfg.radio.sidelobe_floor_dB;

fprintf('\n===== FIG 4.3a - PATRON DE HAZ DEL SATELITE =====\n');
fprintf('perfil            : config_default (Ku de estudio)\n');
fprintf('HPBW              : %.4f deg  (cfg.radio.beamwidth3dB_deg)\n', HPBW);
fprintf('suelo de lobulos  : %.2f dB   (cfg.radio.sidelobe_floor_dB)\n', FLOOR);

%% 1) Muestreo angular: 0 .. 4 anchos de haz (cubre >= 2 lobulos aplanados)
th_max = 4 * HPBW;
th     = linspace(0, th_max, 20001);

G_ideal = beam_gain_dB(th, HPBW);            % sin 3er argumento -> Bessel puro
G_used  = beam_gain_dB(th, HPBW, FLOOR);     % con suelo de lobulos

%% 2) Puntos notables
% (a) media potencia: theta = HPBW/2  ->  -3.0103 dB EXACTOS por construccion
th_3dB = HPBW/2;
G_3dB  = beam_gain_dB(th_3dB, HPBW);

% (b) primer nulo del patron IDEAL: u = 3.8317 (primer cero de 2*J1(u)/u)
u_null1  = 3.8317059702075123;                                  % 1er cero de J1
th_null1 = asind( u_null1 * sind(HPBW/2) / 1.61634 );

% (c) primer angulo en que el SUELO se impone sobre el Bessel
i_floor  = find(G_ideal < FLOOR, 1, 'first');
assert(~isempty(i_floor), 'El suelo nunca se impone en el rango barrido: revisar.');
th_floor = th(i_floor);

%% 3) AUTOVERIFICACION (aborta si algo falla)
err_3dB = abs(G_3dB - (-3.0103));
assert(err_3dB < 1e-3, ...
    'FALLO: beam_gain_dB(HPBW/2,HPBW) = %.6f dB, esperado -3.0103 (err %.2e).', ...
    G_3dB, err_3dB);

min_used = min(G_used);
assert(min_used >= FLOOR - 1e-12, ...
    'FALLO: el patron con suelo baja hasta %.6f dB, por debajo del suelo %.2f dB.', ...
    min_used, FLOOR);

% Coincidencia DENTRO DEL LOBULO PRINCIPAL. Matiz que hay que enunciar bien: el
% lobulo principal llega hasta el primer nulo (th_null1), pero el Bessel ideal ya
% CRUZA el suelo antes (en th_floor) al caer hacia ese nulo, y a partir de ahi las
% dos curvas difieren por definicion. La afirmacion correcta -- y la que se puede
% defender -- es que coinciden EXACTAMENTE mientras el Bessel esta por encima del
% suelo, es decir en 0 <= theta < th_floor, que es la parte util del lobulo
% principal (contiene el pico y los -3 dB).
assert(th_floor < th_null1, ...
    'FALLO: el suelo se impone en %.4f deg, DESPUES del primer nulo %.4f deg.', ...
    th_floor, th_null1);
in_main = th < th_floor;
d_main  = max(abs(G_ideal(in_main) - G_used(in_main)));
assert(d_main < 1e-12, ...
    'FALLO: ideal y con-suelo difieren %.3e dB en el lobulo principal por encima del suelo.', d_main);

fprintf('\n--- autoverificacion ---\n');
fprintf('  G(HPBW/2)                      = %.6f dB   (|err| vs -3.0103 = %.2e)  OK\n', G_3dB, err_3dB);
fprintf('  min del patron con suelo       = %.6f dB   (suelo = %.2f dB)         OK\n', min_used, FLOOR);
fprintf('  max|ideal - suelo| en 0..%.3f deg = %.3e dB                        OK\n', ...
        th_floor, d_main);
fprintf('  (tramo del lobulo principal en que el Bessel esta por encima del suelo)\n');
fprintf('\n--- puntos notables ---\n');
fprintf('  theta_3dB   = HPBW/2      = %.4f deg  ->  %.4f dB\n', th_3dB, G_3dB);
fprintf('  primer nulo (u = %.4f)  = %.4f deg\n', u_null1, th_null1);
fprintf('  el suelo se impone desde   = %.4f deg  (= %.2f anchos de haz)\n', ...
        th_floor, th_floor/HPBW);
fprintf('  rango barrido              = 0 .. %.4f deg (= 4 anchos de haz)\n', th_max);

%% 4) Figura (SIN title; etiquetas ASCII; fuentes >= 11 pt)
fh = figure('Name','Patron de haz del satelite','Color','w','Visible','off');
plot(th, G_ideal, '-',  'LineWidth',1.6, 'Color',[0.00 0.45 0.74]); hold on;
plot(th, G_used,  '--', 'LineWidth',1.8, 'Color',[0.85 0.33 0.10]);
yline(FLOOR, ':', 'Color',[.35 .35 .35], 'LineWidth',1.2, 'HandleVisibility','off');
grid on; box on;
xlim([0 th_max]); ylim([FLOOR-5 0]);
xlabel('Angulo off-boresight [deg]');
ylabel('Ganancia relativa al pico [dB]');

% --- marcas de los tres puntos notables ---
% (1) media potencia. La etiqueta va POR DEBAJO del punto para no chocar con la
%     leyenda, que ocupa la esquina superior derecha.
plot(th_3dB, G_3dB, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6, 'HandleVisibility','off');
text(0.015*th_max, -9.0, ...
     sprintf('HPBW/2 = %.3f deg\n(%.4f dB)', th_3dB, G_3dB), ...
     'FontSize',11, 'VerticalAlignment','top', 'HorizontalAlignment','left');

% (2) primer nulo ideal: la linea vertical marca el angulo; la etiqueta va a su
%     derecha, en la franja vacia entre el lobulo principal y el 1er lobulo lateral.
xline(th_null1, '-.', 'Color',[.35 .35 .35], 'LineWidth',1.0, 'HandleVisibility','off');
text(th_null1+0.072*th_max, -9.0, sprintf('1er nulo ideal:\n%.3f deg', th_null1), ...
     'FontSize',11, 'HorizontalAlignment','left', 'VerticalAlignment','top');

% (3) angulo desde el que el suelo se impone sobre el Bessel
plot(th_floor, FLOOR, 'ks', 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerSize',7, ...
     'HandleVisibility','off');
text(th_floor+0.025*th_max, FLOOR-2.6, ...
     sprintf('el suelo se impone desde %.3f deg', th_floor), ...
     'FontSize',11, 'HorizontalAlignment','left', 'VerticalAlignment','middle');

legend({'patron ideal (Bessel 38.811)', ...
        sprintf('patron usado (suelo %.0f dB)', FLOOR)}, ...
       'Location','northeast', 'FontSize',11);

set(findall(fh,'-property','FontSize'), 'FontSize', 11);
save_fig(fh, OUTDIR, 'fig_4_03a_patron_haz');
close(fh);

fprintf('\n[fig 4.3a] escrita en %s\n', OUTDIR);
