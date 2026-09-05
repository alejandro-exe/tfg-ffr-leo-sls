function users = build_user_grid(cfg)
%BUILD_USER_GRID  Construye la rejilla de usuarios FIJA en tierra (ECEF).
%   Segun cfg.ground.mode:
%     'point'     -> un unico punto (cfg.ground.point = [lat lon])
%     'grid'      -> malla regular en [latlim] x [lonlim] con paso [dlat dlon]
%     'list'      -> lista de puntos (cfg.ground.point = [lat1 lon1; lat2 lon2; ...])
%     'localgrid' -> malla CUADRADA en km alrededor de cfg.ground.point, recortada
%                    a un disco de radio cfg.ground.radius_km con paso
%                    cfg.ground.step_km. Es el modo que necesita la FFR: con un
%                    solo usuario no existen "centro" ni "borde" de celda.
%   Devuelve:
%     users.lat, users.lon  [M x 1]  coordenadas geodesicas
%     users.ecef            [M x 3]  posiciones ECEF [km]
%     users.M                        numero de usuarios
%     users.x_km, users.y_km [M x 1] SOLO en 'localgrid': coordenadas locales
%                    (Este, Norte) en km respecto a cfg.ground.point (para
%                    representar mapas de celdas/colores sin rehacer la conversion)
%
%   AVISO DE COSTE (importante en 'localgrid'): la geometria es un tensor
%   M*N*Nt. Con M ~ cientos de usuarios hay que reducir N (constelacion) y/o Nt
%   (usar cfg.time.dt = 60 s o una ventana corta), o compute_geometry avisara.

switch lower(cfg.ground.mode)
    case 'point'
        lat = cfg.ground.point(1);  lon = cfg.ground.point(2);
    case 'localgrid'
        % Malla local en km (x = Este, y = Norte) centrada en el punto de
        % referencia, recortada al disco de radio radius_km.
        lat0 = cfg.ground.point(1);  lon0 = cfg.ground.point(2);
        Rk   = cfg.ground.radius_km;
        stk  = cfg.ground.step_km;
        nn   = floor(Rk / stk);                 % malla simetrica que incluye el 0
        v    = (-nn:nn) * stk;                  % km
        [X, Y] = meshgrid(v, v);                % X = Este [km], Y = Norte [km]
        keep = (X.^2 + Y.^2) <= Rk^2 + 1e-9;    % recorte circular
        xk = X(keep);  yk = Y(keep);
        % Conversion km -> grados (1 deg de latitud ~ 111.32 km; el meridiano se
        % acorta con el coseno de la latitud del propio punto).
        lat = lat0 + yk / 111.32;
        lon = lon0 + xk ./ (111.32 * cosd(lat));
        users.x_km = xk(:);
        users.y_km = yk(:);
    case 'grid'
        [LON,LAT] = meshgrid( ...
            cfg.ground.lonlim(1):cfg.ground.dlon:cfg.ground.lonlim(2), ...
            cfg.ground.latlim(1):cfg.ground.dlat:cfg.ground.latlim(2));
        lat = LAT(:);  lon = LON(:);
    case 'list'
        lat = cfg.ground.point(:,1);  lon = cfg.ground.point(:,2);
    otherwise
        error('build_user_grid:mode','Modo de terreno desconocido: %s', cfg.ground.mode);
end

users.lat  = lat(:);
users.lon  = lon(:);
users.ecef = geodetic2ecef(cfg, users.lat, users.lon, cfg.ground.alt);
users.M    = numel(users.lat);
end
