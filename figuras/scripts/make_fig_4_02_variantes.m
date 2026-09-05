%% MAKE_FIG_4_02_VARIANTES   TAREA 7 del encargo 2 (+ verificacion de latitud)
%  Genera EN _staging las DOS versiones de la instantanea de constelacion, ambas
%  SIN TITULO, para poder decidir cual va a la memoria:
%
%    fig_4_02_alt_reducida.png  - constelacion REDUCIDA, la que hay ahora en
%                                 figuras/. Sale de geometry_results.mat, cuyo cfg
%                                 trae el "Override de VALIDACION" T=66 / P=6.
%    fig_4_02_alt_completa.png  - PERFIL PRINCIPAL completo, config_default() sin
%                                 tocar: Walker 53:1584/72/1.
%
%  NINGUNA de las dos se promueve: se quedan en _staging para inspeccion.
%  La que ya esta en figuras/ NO se toca.
%
%  REGLA DE RECALCULO: solo re-derivacion GEOMETRICA (build_constellation +
%  propagate + eci2ecef + build_user_grid). No interviene ningun indicador; no se
%  llama a compute_sinr_ffr, compute_interference, compute_kpis, ffr_allocate,
%  ffr_policy, run_one_density ni run_sweep_points.
%
%  ---------------------------------------------------------------------------
%  VERIFICACION DE LA COTA DE LATITUD (seccion 1). Es la pregunta que un tribunal
%  hace mirando esta figura: "por que no hay satelites sobre los polos?". En una
%  orbita CIRCULAR de radio a e inclinacion i, la componente Z esta acotada por
%
%        |Z| <= a * sin(i)        (se alcanza en el punto mas septentrional)
%
%  y por tanto la latitud geocentrica no puede pasar de i. Si la constelacion
%  construida superara esa cota habria un error en build_constellation y esta
%  figura lo estaria enseñando. Se comprueba NUMERICAMENTE, no a ojo: la
%  perspectiva 3D engaña y puntos del hemisferio cercano parecen mas altos.
%  Se comprueba en t0 (lo que la figura dibuja) y sobre una ORBITA COMPLETA.
%  ---------------------------------------------------------------------------
%
%  Uso: run('.../figuras/scripts/make_fig_4_02_variantes.m')

HERE     = fileparts(mfilename('fullpath'));
PROJROOT = fileparts(fileparts(HERE));
addpath(PROJROOT);
STAGE = fullfile(PROJROOT,'REDACCION','_staging');
if ~exist(STAGE,'dir'), mkdir(STAGE); end

fprintf('\n################ TAREA 7: variantes de la figura 4.2 ################\n');

for v = 1:2
    if v == 1
        % --- REDUCIDA: exactamente el cfg con el que se publico la figura ---
        Sg   = load(fullfile(PROJROOT,'geometry_results.mat'), 'cfg','sats','users');
        cfg  = Sg.cfg;   sats = Sg.sats;   users = Sg.users;
        name = 'fig_4_02_alt_reducida';
        etiq = 'REDUCIDA (override de validacion de run_geometry_demo)';
    else
        % --- COMPLETA: perfil principal sin overrides ---
        cfg   = config_default();
        sats  = build_constellation(cfg);
        users = build_user_grid(cfg);
        name  = 'fig_4_02_alt_completa';
        etiq  = 'COMPLETA (perfil principal config_default)';
    end

    K  = cfg.constellations(1);
    t0 = cfg.time.t0;                           % SOLO el instante inicial

    fprintf('\n=================== %s ===================\n', etiq);
    fprintf('  constelacion : %s | T = %d, P = %d, h = %g km, inc = %g deg\n', ...
        K.name, K.T, K.P, K.h, K.inc);

    %% ---------- 1. VERIFICACION DE LA COTA DE LATITUD ----------------------
    a_km  = cfg.const.Re + K.h;                 % orbita circular: a = Re + h
    incd  = K.inc;
    Zmax  = a_km * sind(incd);                  % cota teorica |Z| <= a*sin(i)

    % (a) en t0, que es lo que la figura dibuja
    R0    = propagate(cfg, sats, t0);           % [N x 3 x 1] ECI
    Z0    = R0(:,3,1);
    r0    = sqrt(sum(R0(:,:,1).^2, 2));
    lat0  = asind(Z0 ./ r0);                    % latitud GEOCENTRICA

    % (b) sobre una orbita completa, que es la cota que de verdad importa
    Torb  = 2*pi / sats.n(1);                   % s
    tOrb  = linspace(t0, t0 + Torb, 200);
    Rorb  = propagate(cfg, sats, tOrb);         % [N x 3 x 200]
    Zorb  = squeeze(Rorb(:,3,:));
    rorb  = squeeze(sqrt(sum(Rorb.^2, 2)));
    latOrb= asind(Zorb ./ rorb);

    fprintf('\n  --- cota de latitud ---\n');
    fprintf('    a = Re + h                 = %.3f km\n', a_km);
    fprintf('    cota teorica |Z| = a*sin(i)= %.3f km\n', Zmax);
    fprintf('    max|Z| medido en t0        = %.3f km   (margen %+.3f km)\n', ...
            max(abs(Z0)), Zmax - max(abs(Z0)));
    fprintf('    max|Z| sobre 1 orbita      = %.3f km   (margen %+.6f km)\n', ...
            max(abs(Zorb(:))), Zmax - max(abs(Zorb(:))));
    fprintf('    max|lat| geocentrica en t0 = %.4f deg  (inclinacion = %.4f deg)\n', ...
            max(abs(lat0)), incd);
    fprintf('    max|lat| sobre 1 orbita    = %.4f deg  (margen %+.2e deg)\n', ...
            max(abs(latOrb(:))), incd - max(abs(latOrb(:))));
    fprintf('    radio orbital: min %.4f / max %.4f km (circular: dispersion %.2e km)\n', ...
            min(rorb(:)), max(rorb(:)), max(rorb(:)) - min(rorb(:)));

    tolZ = 1e-6;                                 % 1 mm: solo precision de maquina
    assert(max(abs(Zorb(:))) <= Zmax + tolZ, ...
        ['FALLO: |Z| llega a %.3f km y la cota de una orbita circular de %.3f km ' ...
         'con i = %g deg es %.3f km. Habria un error en build_constellation.'], ...
        max(abs(Zorb(:))), a_km, incd, Zmax);
    assert(max(abs(latOrb(:))) <= incd + 1e-9, ...
        'FALLO: la latitud geocentrica alcanza %.6f deg > inclinacion %.6f deg.', ...
        max(abs(latOrb(:))), incd);
    fprintf('    VERIFICACION OK: la cota se respeta, el clareo polar es REAL.\n');

    %% ---------- 2. Indice de plano orbital (para colorear uno) -------------
    % Se agrupa por RAAN en vez de asumir el orden de build_constellation: es
    % robusto aunque cambie el recorrido de los bucles.
    [~, ~, planeIdx] = unique(round(sats.raan, 9));
    nPl = max(planeIdx);
    fprintf('\n  --- estructura de planos ---\n');
    fprintf('    planos detectados por RAAN = %d (esperado P = %d)\n', nPl, K.P);
    assert(nPl == K.P, 'El numero de RAAN distintos (%d) no coincide con P (%d).', nPl, K.P);
    selP = 1;                                    % plano que se resalta
    inP  = (planeIdx == selP);
    fprintf('    satelites por plano        = %d (esperado T/P = %d)\n', nnz(inP), K.T/K.P);
    fprintf('    plano resaltado            = %d de %d, RAAN = %.2f deg\n', ...
            selP, nPl, rad2deg(sats.raan(find(inP,1))));

    %% ---------- 3. Figura --------------------------------------------------
    R_ecef = eci2ecef(cfg, R0, t0);
    P0     = R_ecef(:,:,1);

    f = figure('Name','Constelacion 3D','Color','w','Visible','off');
    hold on; axis equal; grid on;

    % Esfera terrestre SIN brillo especular: relleno plano y uniforme. El
    % resplandor de camlight+gouraud no aporta forma y en escala de grises se
    % come los satelites que caen encima.
    [xe,ye,ze] = sphere(48); Re = cfg.const.Re;
    surf(Re*xe, Re*ye, Re*ze, 'FaceColor',[0.86 0.89 0.94], 'EdgeColor','none', ...
         'FaceAlpha',0.55, 'FaceLighting','none');

    % Satelites: el resto de la constelacion primero, el plano resaltado encima.
    msz = 12;  if v == 2, msz = 5; end
    scatter3(P0(~inP,1), P0(~inP,2), P0(~inP,3), msz, ...
             'MarkerFaceColor',[0.30 0.55 0.80], 'MarkerEdgeColor','none', ...
             'MarkerFaceAlpha',0.65);
    scatter3(P0(inP,1),  P0(inP,2),  P0(inP,3),  msz*2.6, ...
             'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','none');

    % Punto de observacion: con 1584 puntos azules, la estrella se perdia.
    gu = users.ecef(1,:);
    scatter3(gu(1),gu(2),gu(3), 320, 'p', 'MarkerFaceColor',[0.90 0.10 0.10], ...
             'MarkerEdgeColor','k', 'LineWidth',1.3);

    xlabel('X_{ECEF} [km]'); ylabel('Y_{ECEF} [km]'); zlabel('Z_{ECEF} [km]');
    % sin title: el pie va en LaTeX
    % Leyenda FUERA del eje: al noreste tapaba justo el casquete polar norte, que
    % es precisamente el rasgo que la figura tiene que enseñar (el clareo polar).
    legend({'esfera terrestre', ...
            sprintf('resto de la constelacion (%d de %d)', nnz(~inP), K.T), ...
            sprintf('un plano orbital (%d satelites)', nnz(inP)), ...
            'punto de observacion'}, ...
           'Location','southoutside', 'NumColumns',2);
    view(35,20);
    set(findall(f,'-property','FontSize'),'FontSize',11);
    save_fig(f, STAGE, name);
    close(f);
    fprintf('\n  -> %s.png (en _staging, NO promovida)\n', name);
end

fprintf('\n[tarea 7] las dos variantes quedan en %s para decidir.\n', STAGE);
