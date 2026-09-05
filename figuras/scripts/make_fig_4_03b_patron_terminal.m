%% MAKE_FIG_4_03B_PATRON_TERMINAL   Candidata B a la figura 4.3
%  Envolvente de ganancia del TERMINAL en tierra (ITU-R S.1428-1), en dBi
%  ABSOLUTOS, frente al angulo respecto al eje de punteria, de 0 a 180 deg.
%  Eje X LOGARITMICO: es la representacion habitual de S.1428 y la unica forma de
%  ver a la vez el lobulo principal (decimas de grado) y la cola (decenas).
%
%  PERFIL: config_default() (perfil Ku de ESTUDIO). NUNCA config_calib.
%
%  Que ilustra: la discriminacion del VSAT fuera de eje (35-48 dB), que es la
%  MITAD de la "doble discriminacion" con la que el TFG explica por que la
%  interferencia INTER-satelite resulta despreciable.
%
%  AVISO ESPERADO Y DELIBERADAMENTE NO SILENCIADO: con este perfil
%  D/lambda ~ 23.9 cae en el tramo 20 <= D/L <= 25 de S.1428-1 mientras el codigo
%  implementa el de 25 < D/L <= 100. term_gain_dB emite warning('term_gain_dB:regime25').
%  Es CONSERVADOR (sobreestima 5 dB la cola por encima de 120 deg). Se captura y
%  se reporta, no se corrige.
%
%  SIN TITULO (el pie va en LaTeX). Solo lee el motor; escribe unicamente en
%  REDACCION/figuras/.
%
%  Uso: desde la raiz del proyecto,
%       run('figuras/scripts/make_fig_4_03b_patron_terminal.m')

% Rutas ancladas a la ubicacion del script (ver nota en regen_fig_4_02_*).
HERE     = fileparts(mfilename('fullpath'));
PROJROOT = fileparts(fileparts(HERE));
addpath(PROJROOT);
OUTDIR   = fullfile(PROJROOT,'REDACCION','figuras');

cfg  = config_default();
Gmax = cfg.radio.Grx_max_dBi;
eta  = cfg.radio.term_efficiency;
fGHz = cfg.radio.freq_GHz;

%% 0) Rearmar el aviso de regimen (term_gain_dB avisa UNA vez por (Gmax,eta))
clear term_gain_dB
lastwarn('');

%% 1) Parametros derivados de la envolvente (mismas formulas que term_gain_dB)
DL    = sqrt( 10^(Gmax/10) / (eta * pi^2) );   % D/lambda derivado de Gmax y eta
G1    = 29 - 25*log10(95/DL);                  % meseta
phi_m = (20/DL) * sqrt(max(Gmax - G1, 0));     % fin del lobulo principal
phi_r = 95/DL;                                 % inicio del tramo 29-25log10(phi)

lambda_m = 299792458 / (fGHz*1e9);             % m
D_geom_m = DL * lambda_m;                      % diametro implicado por D/lambda

fprintf('\n===== FIG 4.3b - PATRON DEL TERMINAL (ITU-R S.1428-1) =====\n');
fprintf('perfil                 : config_default (Ku de estudio)\n');
fprintf('--- entradas ---\n');
fprintf('  Grx_max              = %.4f dBi   (cfg.radio.Grx_max_dBi)\n', Gmax);
fprintf('  eficiencia apertura  = %.4f       (cfg.radio.term_efficiency)\n', eta);
fprintf('  frecuencia           = %.4f GHz   -> lambda = %.6f m\n', fGHz, lambda_m);
fprintf('--- parametros derivados de la envolvente ---\n');
fprintf('  D/lambda             = %.4f\n', DL);
fprintf('  (diametro implicado) = %.4f m\n', D_geom_m);
fprintf('  G1 (meseta)          = %.4f dBi\n', G1);
fprintf('  phi_m (fin lob.ppal) = %.4f deg\n', phi_m);
fprintf('  phi_r (inicio 29-25log) = %.4f deg\n', phi_r);

%% 2) AUTOVERIFICACION: en phi = 0 la ganancia debe ser EXACTAMENTE Grx_max
G0 = term_gain_dB(cfg, 0);
assert(G0 == Gmax, ...
    'FALLO: term_gain_dB(cfg,0) = %.12f dBi, esperado exactamente %.12f dBi.', G0, Gmax);
fprintf('--- autoverificacion ---\n');
fprintf('  term_gain_dB(cfg,0)  = %.12f dBi  ==  Grx_max  OK (igualdad exacta)\n', G0);

%% 3) Muestreo en escala logaritmica (el 0 no cabe en un eje log: se arranca en 0.01)
phi = logspace(log10(0.01), log10(180), 20001);
G   = term_gain_dB(cfg, phi);

[wmsg, wid] = lastwarn;
if ~isempty(wmsg)
    fprintf('--- aviso de MATLAB capturado (esperado, NO silenciado) ---\n');
    fprintf('  identificador: %s\n', wid);
    fprintf('  texto        : %s\n', wmsg);
end

%% 4) Figura (SIN title; etiquetas ASCII; fuentes >= 11 pt)
fh = figure('Name','Patron del terminal S.1428-1','Color','w','Visible','off');
semilogx(phi, G, '-', 'LineWidth',1.8, 'Color',[0.00 0.45 0.74]); hold on;
grid on; box on;
xlim([0.01 200]); ylim([-15 Gmax+5]);   % 200 deg de tope: margen para las etiquetas
xlabel('Angulo respecto al eje de punteria [deg]');
ylabel('Ganancia del terminal [dBi]');

% --- transiciones de la envolvente ---
%  phi_m y phi_r caen casi encima (3.856 vs 3.971 deg): en eje log sus etiquetas
%  se solapan, asi que van JUNTAS en un bloque unico arriba a la derecha. Los dos
%  escalones de cola llevan etiqueta VERTICAL sobre su propia linea, que es lo
%  unico que cabe sin salirse del lienzo.
trans = [phi_m, phi_r, 33.1, 80];
for k = 1:numel(trans)
    xline(trans(k), '--', 'Color',[.45 .45 .45], 'LineWidth',1.1, 'HandleVisibility','off');
    plot(trans(k), term_gain_dB(cfg, trans(k)), 'ko', 'MarkerFaceColor','k', ...
         'MarkerSize',5, 'HandleVisibility','off');
end

% Bloque phi_m / phi_r: va a la IZQUIERDA, en la franja vacia bajo la meseta de
% 35.3 dBi (la curva es plana hasta ~2 deg, asi que ahi no hay nada que tapar).
text(0.012, 31, ...
     sprintf('phi_m = %.3f deg : fin del lobulo principal\nphi_r = %.3f deg : inicio de 29-25log(phi)', ...
             phi_m, phi_r), ...
     'FontSize',11, 'VerticalAlignment','top', 'HorizontalAlignment','left');

% Escalones de cola: etiqueta VERTICAL que arranca justo encima de su tramo y
% crece hacia arriba, por la zona vacia que queda sobre la cola.
text(36, -1, '33.1 deg : escalon a -9 dBi', 'FontSize',11, ...
     'Rotation',90, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
text(87, -1, '80 deg : escalon a -4 dBi', 'FontSize',11, ...
     'Rotation',90, 'HorizontalAlignment','left', 'VerticalAlignment','middle');

% meseta G1
yline(G1, ':', 'Color',[0.85 0.33 0.10], 'LineWidth',1.2, 'HandleVisibility','off');
text(0.012, G1+1.6, sprintf('G1 = %.2f dBi', G1), 'FontSize',11, ...
     'Color',[0.85 0.33 0.10]);

legend({sprintf('S.1428-1, D/lambda = %.2f', DL)}, 'Location','southwest', 'FontSize',11);

set(findall(fh,'-property','FontSize'), 'FontSize', 11);
save_fig(fh, OUTDIR, 'fig_4_03b_patron_terminal');
close(fh);

fprintf('\n[fig 4.3b] escrita en %s\n', OUTDIR);
