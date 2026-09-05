%% REGEN_FIG_4_02_CONSTELACION_ECEF
%  Regenera la figura 4.2 (instantanea 3D de la constelacion en t0, marco ECEF)
%  SIN TITULO, para que el pie vaya en LaTeX.
%
%  ORIGEN: run_geometry_demo.m, seccion "9) Figura 6", lineas 108-117.
%  El graficado se reproduce LINEA A LINEA; lo unico que se omite es title(...).
%
%  NOTA IMPORTANTE (declarada en INFORME_FIGURAS.md):
%  geometry_results.mat guarda {cfg, tvec, sats, users, G, S} pero NO guarda
%  R_ecef, que es lo que se dibuja. Por eso este script re-deriva las posiciones
%  con propagate + eci2ecef A PARTIR DEL cfg Y sats GUARDADOS. No es re-lanzar el
%  experimento: es cinematica kepleriana determinista sobre la misma entrada, y no
%  toca ninguna funcion de fisica ni de KPI (no llama a compute_sinr_ffr,
%  compute_interference, compute_kpis, run_one_density ni run_sweep_points).
%  La fidelidad se comprueba a posteriori comparando pixel a pixel la variante
%  ESTRICTA contra figs_geom/geom_f_constelacion_3d_ecef_t0.png.
%
%  SALIDAS (ambas dentro de REDACCION/, no escribe en ninguna otra ruta):
%    REDACCION/figuras/fig_4_02_constelacion_ecef.png       <- entregable (fuentes >= 11 pt)
%    REDACCION/figuras/_verificacion/fig_4_02_estricta.png  <- control de fidelidad
%
%  Uso: desde la raiz del proyecto,
%       run('figuras/scripts/regen_fig_4_02_constelacion_ecef.m')

% Rutas ANCLADAS A LA UBICACION DEL SCRIPT, no al directorio de trabajo: `run`
% cambia el cwd a la carpeta del fichero, y con rutas relativas las salidas
% acababan en una carpeta anidada bajo el propio directorio del script.
HERE     = fileparts(mfilename('fullpath'));       % .../figuras/scripts
PROJROOT = fileparts(fileparts(HERE));             % raiz del proyecto
addpath(PROJROOT);
OUTDIR   = fullfile(PROJROOT,'REDACCION','figuras');
VERDIR   = fullfile(OUTDIR,'_verificacion');

%% 1) Cargar el .mat publicado y re-derivar SOLO la geometria orbital
Sg   = load(fullfile(PROJROOT,'geometry_results.mat'), 'cfg', 'tvec', 'sats', 'users', 'G');
cfg  = Sg.cfg;  sats = Sg.sats;  users = Sg.users;

fprintf('[fig 4.2] cfg del .mat: T=%d, P=%d, N=%d satelites, Nt=%d\n', ...
    cfg.constellations(1).T, cfg.constellations(1).P, sats.N, numel(Sg.tvec));

tvec = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
assert(isequal(size(tvec), size(Sg.tvec)) && max(abs(tvec - Sg.tvec)) == 0, ...
    'El tvec reconstruido no coincide con el guardado: PARAR y revisar.');

R_eci  = propagate(cfg, sats, tvec);
R_ecef = eci2ecef(cfg, R_eci, tvec);

%% 1b) PRUEBA DE QUE EL R_ecef RE-DERIVADO ES EL MISMO QUE EL PUBLICADO
%  La comparacion pixel a pixel de la figura no puede distinguir una diferencia de
%  DATOS de una de RASTERIZADO. Esta si: geometry_results.mat SI guarda G
%  (elev/azimut/rango/visibilidad del usuario a cada satelite en cada instante),
%  que es funcion directa de R_ecef. Si G recalculado coincide BIT A BIT con el
%  guardado, el R_ecef re-derivado es exactamente el que produjo la figura.
%  compute_geometry es geometria pura (no es fisica ni KPI) y aqui se usa SOLO
%  como verificacion: no interviene en el graficado.
Gchk = compute_geometry(cfg, users, R_ecef);
campos = {'el','az','range','vis'};
fprintf('[fig 4.2] verificacion de datos (G recalculado vs G guardado):\n');
for f = 1:numel(campos)
    c = campos{f};
    if ~isfield(Sg.G, c) || ~isfield(Gchk, c), continue; end
    d = max(abs(double(Gchk.(c)(:)) - double(Sg.G.(c)(:))));
    fprintf('    max|dif| en G.%-6s = %.3e\n', c, d);
    assert(d == 0, ...
        'FALLO: G.%s recalculado difiere del guardado (max|dif| = %.3e). PARAR.', c, d);
end
fprintf('    -> R_ecef re-derivado IDENTICO al publicado (max|dif| = 0 en todos los campos).\n');

%% 2) Graficado (copia literal de run_geometry_demo.m:109-116)
%  Se dibuja TRES veces:
%    v=1 REPLICA  - identica al original, TITULO INCLUIDO. Sirve para el control de
%                   fidelidad: se compara pixel a pixel contra el PNG publicado. Es
%                   la unica comparacion limpia, porque quitar el titulo cambia el
%                   recorte de exportgraphics y desplaza/estrecha el lienzo.
%    v=2 ESTRICTA - igual que la replica pero SIN title (fuentes por defecto).
%    v=3 ENTREGABLE - sin title y con fuentes >= 11 pt (regla de estilo).
for v = 1:3
    fh = figure('Name','Constelacion 3D','Color','w','Visible','off');
    hold on; axis equal; grid on;
    [xe,ye,ze] = sphere(36); Re = cfg.const.Re;
    surf(Re*xe, Re*ye, Re*ze, 'FaceColor',[.85 .9 1], 'EdgeColor','none','FaceAlpha',.6);
    P0 = R_ecef(:,:,1);
    scatter3(P0(:,1),P0(:,2),P0(:,3),12,'filled');
    gu = users.ecef(1,:); scatter3(gu(1),gu(2),gu(3),80,'rp','filled');
    xlabel('X_{ECEF} [km]'); ylabel('Y_{ECEF} [km]'); zlabel('Z_{ECEF} [km]');
    if v == 1
        title('Constelacion en t_0 (ECEF)');   % SOLO en la replica de control
    end
    % en v=2 y v=3 el title se OMITE a proposito: el pie va en LaTeX.
    view(35,20); camlight; lighting gouraud;

    switch v
        case 1, save_fig(fh, VERDIR, 'fig_4_02_replica');
        case 2, save_fig(fh, VERDIR, 'fig_4_02_estricta');
        case 3
            set(findall(fh,'-property','FontSize'), 'FontSize', 11);
            save_fig(fh, OUTDIR, 'fig_4_02_constelacion_ecef');
    end
    close(fh);
end

fprintf('[fig 4.2] escritas en %s y %s\n', OUTDIR, VERDIR);
