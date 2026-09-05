%% REGEN_FIG_4_05_CELDAS_COLOREADO
%  Regenera la figura 4.5 (mapa de celdas Earth-fixed, coloreado Delta=3 y rejilla
%  de usuarios) SIN TITULO, para que el pie vaya en LaTeX.
%
%  ORIGEN: run_ffr_demo.m, seccion "(e)", lineas 565-580.
%  El graficado se reproduce LINEA A LINEA; lo unico que se omite es title(...).
%
%  NOTA IMPORTANTE (declarada en INFORME_FIGURAS.md):
%  ffr_results.mat guarda BL (el layout de haces, que es lo que fija las celdas)
%  pero NO guarda `users`. Por eso este script reconstruye la rejilla con
%  build_user_grid(cfg) usando el cfg GUARDADO: build_user_grid es un generador
%  determinista de malla a partir de cfg.ground (mode/radius_km/step_km), sin
%  fisica ni aleatoriedad. Igual que ffr_coloring(BL,3), que es combinatoria pura
%  sobre los indices axiales de BL. NO se llama a compute_sinr_ffr,
%  compute_interference, compute_kpis, run_one_density ni run_sweep_points.
%
%  SALIDAS (ambas dentro de REDACCION/, no escribe en ninguna otra ruta):
%    REDACCION/figuras/fig_4_05_celdas_coloreado.png        <- entregable (fuentes >= 11 pt)
%    REDACCION/figuras/_verificacion/fig_4_05_estricta.png  <- control de fidelidad
%
%  Uso: desde la raiz del proyecto,
%       run('figuras/scripts/regen_fig_4_05_celdas_coloreado.m')

% Rutas ANCLADAS A LA UBICACION DEL SCRIPT, no al directorio de trabajo: `run`
% cambia el cwd a la carpeta del fichero, y con rutas relativas las salidas
% acababan en una carpeta anidada bajo el propio directorio del script.
HERE     = fileparts(mfilename('fullpath'));       % .../figuras/scripts
PROJROOT = fileparts(fileparts(HERE));             % raiz del proyecto
addpath(PROJROOT);
OUTDIR   = fullfile(PROJROOT,'REDACCION','figuras');
VERDIR   = fullfile(OUTDIR,'_verificacion');

%% 1) Cargar el .mat publicado y reconstruir la rejilla de usuarios
Sf  = load(fullfile(PROJROOT,'ffr_results.mat'), 'cfg', 'BL');
cfg = Sf.cfg;  BL = Sf.BL;

fprintf(['[fig 4.5] cfg del .mat: ground.mode=%s, radius=%g km, step=%g km | ' ...
         'BL: nBeams=%d, s=%.4f km, lattice=%s\n'], ...
    cfg.ground.mode, cfg.ground.radius_km, cfg.ground.step_km, ...
    BL.nBeams, BL.spacing_km, BL.lattice);

users = build_user_grid(cfg);
fprintf('[fig 4.5] rejilla reconstruida: M = %d usuarios\n', users.M);

colD = ffr_coloring(BL, 3);

%% 2) Graficado (copia literal de run_ffr_demo.m:566-579)
%  Se dibuja TRES veces:
%    v=1 REPLICA  - identica al original, TITULO INCLUIDO, para el control de
%                   fidelidad contra el PNG publicado. Hace falta porque el titulo
%                   de esta figura es MAS ANCHO que los ejes, asi que al quitarlo
%                   exportgraphics recorta un lienzo mas estrecho y la comparacion
%                   directa sin titulo seria imposible.
%    v=2 ESTRICTA - sin title, fuentes por defecto.
%    v=3 ENTREGABLE - sin title y con fuentes >= 11 pt (regla de estilo).
for v = 1:3
    fh = figure('Name','Celdas Earth-fixed y coloreado FFR','Color','w','Visible','off');
    cmap = lines(3);
    th = linspace(0,2*pi,40);
    Rc = BL.spacing_km/sqrt(3);                       % circunradio de la celda hex
    for b = 1:BL.nBeams
        p = BL.cell_km(b,:);
        fill(p(1)+Rc*cos(th), p(2)+Rc*sin(th), cmap(colD(b),:), ...
            'FaceAlpha',0.25, 'EdgeColor',[.6 .6 .6]); hold on;
    end
    plot(users.x_km, users.y_km, 'k.', 'MarkerSize',4);
    axis equal; grid on; xlabel('Este [km]'); ylabel('Norte [km]');
    if v == 1
        title(sprintf(['Celdas Earth-fixed (%d haces, s=%.1f km), coloreado ' ...
                       '\\Delta=3 y rejilla de usuarios'], BL.nBeams, BL.spacing_km));
    end
    % en v=2 y v=3 el title se OMITE a proposito: el pie va en LaTeX.

    switch v
        case 1, save_fig(fh, VERDIR, 'fig_4_05_replica');
        case 2, save_fig(fh, VERDIR, 'fig_4_05_estricta');
        case 3
            set(findall(fh,'-property','FontSize'), 'FontSize', 11);
            save_fig(fh, OUTDIR, 'fig_4_05_celdas_coloreado');
    end
    close(fh);
end

fprintf('[fig 4.5] escritas en %s y %s\n', OUTDIR, VERDIR);
