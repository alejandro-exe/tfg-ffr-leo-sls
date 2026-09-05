%% RUN_GEOMETRY_DEMO  Demostrador y validacion de la capa geometrica del SLS LEO-FFR.
%  Flujo:  config -> constelacion -> propagacion (ECI) -> ECEF -> usuarios ->
%          geometria (elev/az/rango/visibilidad) -> satelite servidor -> salidas.
%
%  Esta capa NO calcula aun SINR ni FFR: produce las ENTRADAS que consumiran
%  los modulos de radioenlace (rango, elevacion -> FSPL y atmosfera), de
%  interferencia (conjunto de satelites visibles co-canal) y de FFR
%  (clasificacion centro/borde a partir de la geometria del servidor).
%
%  Uso:  abrir esta carpeta en MATLAB y ejecutar.  Requiere solo MATLAB base.

clear; clc; close all;

FIGDIR = 'figs_geom';        % PNG a 300 dpi via save_fig (fuente unica)

cfg = config_default();

% --- Override de VALIDACION: subconjunto pequeno para depurar y graficar ---
%     (la constelacion real 1584/72 da trazas ilegibles y es lenta de depurar;
%      comenta este bloque para usar los valores reales de config_default).
cfg.constellations(1).T = 66;     % 6 planos x 11 satelites
cfg.constellations(1).P = 6;
% ---------------------------------------------------------------------------

%% 1) Construir y propagar
tvec  = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;   % 1 x Nt
sats  = build_constellation(cfg);
R_eci = propagate(cfg, sats, tvec);
R_ecef= eci2ecef(cfg, R_eci, tvec);
users = build_user_grid(cfg);

%% 2) Geometria y satelite servidor
G = compute_geometry(cfg, users, R_ecef);
S = associate_serving(cfg, G);

%% 3) Comprobaciones por consola (validacion teorica)
fprintf('\n===== VALIDACION DE LA GEOMETRIA =====\n');
fprintf('Satelites totales: %d   |  Instantes: %d  (dt=%g s, %g min)\n', ...
        sats.N, numel(tvec), cfg.time.dt, cfg.time.duration/60);
for c = 1:numel(cfg.constellations)
    K = cfg.constellations(c);
    a = cfg.const.Re + K.h;  n = sqrt(cfg.const.mu/a^3);  Tmin = 2*pi/n/60;
    fprintf('  [%s] h=%g km  inc=%g deg  ->  periodo teorico = %.2f min\n', ...
            K.name, K.h, K.inc, Tmin);
end
fprintf('Usuario 1 (lat=%.2f, lon=%.2f), elev. minima = %g deg:\n', ...
        users.lat(1), users.lon(1), cfg.geom.minElev);
fprintf('   satelites visibles: min=%d  medio=%.2f  max=%d\n', ...
        min(S.nVis(1,:)), mean(S.nVis(1,:)), max(S.nVis(1,:)));
fprintf('   cobertura temporal: %.1f%% de los instantes con >=1 satelite\n', ...
        100*mean(S.nVis(1,:) >= 1));
handovers = sum(diff(S.idx(1,:)) ~= 0 & ~isnan(diff(S.idx(1,:))));
fprintf('   cambios de satelite servidor (handovers): %d\n', handovers);
fprintf('======================================\n\n');

%% 4) Figura 1 - Traza terrestre (sub-satellite points)
rmag = sqrt(sum(R_ecef.^2,2));                         % N x 1 x Nt
subLat = squeeze(asind( R_ecef(:,3,:) ./ rmag ));      % N x Nt
subLon = squeeze(atan2d( R_ecef(:,2,:), R_ecef(:,1,:)));
f1 = figure('Name','Traza terrestre','Color','w');
plot(subLon(:), subLat(:), '.', 'MarkerSize', 3); hold on;
plot(users.lon, users.lat, 'rp', 'MarkerFaceColor','r', 'MarkerSize',12);
xlabel('Longitud [deg]'); ylabel('Latitud [deg]');
xlim([-180 180]); ylim([-90 90]); grid on;
title(sprintf('Traza terrestre (|lat| acotada por inc=%g deg)', cfg.constellations(1).inc));
legend('Sub-satellite points','Usuario','Location','southoutside','Orientation','horizontal');
save_fig(f1, FIGDIR, 'geom_a_traza_terrestre_walker');

%% 5) Figura 2 - Satelites visibles vs tiempo (usuario 1)
f2 = figure('Name','Satelites visibles','Color','w');
stairs(tvec/60, S.nVis(1,:), 'LineWidth',1.2); grid on;
xlabel('Tiempo [min]'); ylabel('N.o de satelites visibles');
title('Satelites visibles sobre el usuario (futuros interferentes co-canal)');
save_fig(f2, FIGDIR, 'geom_b_satelites_visibles_vs_tiempo');

%% 6) Figura 3 - Elevacion del satelite servidor vs tiempo (vista "desde el suelo")
f3 = figure('Name','Elevacion del servidor','Color','w');
plot(tvec/60, S.el(1,:), 'LineWidth',1.4); grid on;
yline(cfg.geom.minElev,'--','Elev. minima');
xlabel('Tiempo [min]'); ylabel('Elevacion del servidor [deg]'); ylim([0 90]);
title('Arco de elevacion del satelite servidor (sube-pico-baja; saltos = handover)');
save_fig(f3, FIGDIR, 'geom_c_elevacion_servidor_vs_tiempo');

%% 7) Figura 4 - Rango de inclinacion del servidor vs tiempo
f4 = figure('Name','Rango del servidor','Color','w');
plot(tvec/60, S.range(1,:), 'LineWidth',1.4); grid on;
xlabel('Tiempo [min]'); ylabel('Slant range [km]');
title('Distancia al satelite servidor (minima en el cenit; fija la FSPL)');
save_fig(f4, FIGDIR, 'geom_d_rango_servidor_vs_tiempo');

%% 8) Figura 5 - Skyplot (azimut/elevacion) de los pases sobre el usuario 1
elU = squeeze(G.el(1,:,:));  azU = squeeze(G.az(1,:,:));  viU = squeeze(G.vis(1,:,:));
ev = elU(viU);  av = azU(viU);
rr = 90 - ev;                       % radio del skyplot: cenit en el centro
xx = rr.*sind(av);  yy = rr.*cosd(av);
f5 = figure('Name','Skyplot','Color','w'); hold on; axis equal; axis off;
th = linspace(0,2*pi,200);
for e = [0 30 60]                   % circulos de elevacion
    plot((90-e)*cos(th),(90-e)*sin(th),'Color',[.8 .8 .8]);
    text(0,(90-e),sprintf(' %d',e),'Color',[.5 .5 .5],'FontSize',8);
end
plot([0 0],[-90 90],'Color',[.85 .85 .85]); plot([-90 90],[0 0],'Color',[.85 .85 .85]);
scatter(xx, yy, 8, 'filled');
text(0,95,'N','HorizontalAlignment','center'); text(95,0,'E','HorizontalAlignment','center');
title('Skyplot: pases visibles sobre el usuario (centro = cenit)');
save_fig(f5, FIGDIR, 'geom_e_skyplot_pases_visibles');

%% 9) Figura 6 - Instantanea 3D de la constelacion en t0
f6 = figure('Name','Constelacion 3D','Color','w'); hold on; axis equal; grid on;
[xe,ye,ze] = sphere(36); Re = cfg.const.Re;
surf(Re*xe, Re*ye, Re*ze, 'FaceColor',[.85 .9 1], 'EdgeColor','none','FaceAlpha',.6);
P0 = R_ecef(:,:,1);
scatter3(P0(:,1),P0(:,2),P0(:,3),12,'filled');
gu = users.ecef(1,:); scatter3(gu(1),gu(2),gu(3),80,'rp','filled');
xlabel('X_{ECEF} [km]'); ylabel('Y_{ECEF} [km]'); zlabel('Z_{ECEF} [km]');
title('Constelacion en t_0 (ECEF)'); view(35,20); camlight; lighting gouraud;
save_fig(f6, FIGDIR, 'geom_f_constelacion_3d_ecef_t0');

fprintf('\n[figuras] 6 PNG a 300 dpi en %s%s\n', FIGDIR, filesep);

%% 10) Guardar resultados para los modulos posteriores
save('geometry_results.mat','cfg','tvec','sats','users','G','S');
fprintf('Resultados guardados en geometry_results.mat\n');
