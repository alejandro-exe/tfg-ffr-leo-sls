function G = compute_geometry(cfg, users, R_ecef, minElev_sat)
%COMPUTE_GEOMETRY  Geometria usuario-satelite desde el suelo (topocentrica ENU).
%   G = compute_geometry(cfg, users, R_ecef)
%   G = compute_geometry(cfg, users, R_ecef, minElev_sat)
%
%   Para cada usuario m, satelite j e instante k calcula:
%     G.el    [M x N x Nt]  elevacion  [deg]
%     G.az    [M x N x Nt]  azimut     [deg, desde el Norte hacia el Este]
%     G.range [M x N x Nt]  rango de inclinacion (slant range) [km]
%     G.vis   [M x N x Nt]  logico: elevacion >= mascara del satelite j
%
%   MASCARA POR SATELITE (minElev_sat, OPCIONAL [N x 1] deg): elevacion minima de
%   OPERACION de cada satelite, tal como la devuelve build_constellation en
%   sats.minElev. Si se omite (o va vacia) se usa el escalar cfg.geom.minElev para
%   TODOS, que es el comportamiento historico -- los llamadores de una sola
%   constelacion quedan bit-identicos.
%
%   POR QUE HACE FALTA: los operadores declaran elevaciones minimas DISTINTAS
%   (Starlink 25 deg, OneWeb 55 deg). Con una unica mascara global, en un escenario
%   multi-operador entrarian como interferentes satelites ajenos a elevaciones a las
%   que NO transmiten, sobreestimando la probabilidad de evento in-line de E4. La
%   mascara filtra a la vez la seleccion de servidor (associate_serving) y el conjunto
%   de interferentes (compute_interference), que es justo lo que se quiere.
%
%   NOTA: cfg.geom.minElev sigue siendo el suelo de la tabla P.618 de atm_loss_dB y su
%   clave de cache. Debe ser <= min(minElev_sat) para que no haya que extrapolar.
%
%   ENU: se proyecta el vector (sat - usuario) en el sistema local Este-Norte-Up
%   del usuario; elevacion = asin(Up/range), azimut = atan2(Este, Norte).
%
%   AVISO de memoria: los arrays son M*N*Nt. Para validar usa 'point' (M=1) y
%   un subconjunto de satelites; sube a la constelacion completa despues.

M = users.M;  [N,~,Nt] = size(R_ecef);

if M*N*Nt > 5e7
    warning('compute_geometry:size', ...
        'Tensor grande (%d elementos). Reduce usuarios/satelites o trocea en tiempo.', M*N*Nt);
end

G.el    = zeros(M,N,Nt);
G.az    = zeros(M,N,Nt);
G.range = zeros(M,N,Nt);

lat = deg2rad(users.lat);
lon = deg2rad(users.lon);

for m = 1:M
    g  = users.ecef(m,:);                    % 1 x 3 (ECEF del usuario)
    sl = sin(lat(m)); cl = cos(lat(m));
    so = sin(lon(m)); co = cos(lon(m));
    Renu = [ -so      ,  co      , 0 ; ...    % fila Este
             -sl*co   , -sl*so   , cl; ...    % fila Norte
              cl*co   ,  cl*so   , sl ];       % fila Up

    for k = 1:Nt
        d   = R_ecef(:,:,k) - g;             % N x 3  (vector usuario -> satelite)
        enu = (Renu * d.').';                % N x 3  en coordenadas ENU
        rng = sqrt(sum(enu.^2, 2));
        G.range(m,:,k) = rng;
        G.el(m,:,k)    = asind( enu(:,3) ./ rng );
        G.az(m,:,k)    = mod( atan2d(enu(:,1), enu(:,2)), 360 );
    end
end

% Mascara de visibilidad: por satelite si se ha dado, global si no.
if nargin < 4 || isempty(minElev_sat)
    G.vis = G.el >= cfg.geom.minElev;
else
    if numel(minElev_sat) ~= N
        error('compute_geometry:minElevSize', ...
            'minElev_sat tiene %d elementos y hay %d satelites.', numel(minElev_sat), N);
    end
    % [M x N x Nt] >= [1 x N] por expansion implicita en la dimension de satelites
    G.vis = G.el >= reshape(minElev_sat, 1, N);
end
end
