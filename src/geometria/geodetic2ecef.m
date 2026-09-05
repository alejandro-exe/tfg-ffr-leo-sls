function r = geodetic2ecef(cfg, lat_deg, lon_deg, h_km)
%GEODETIC2ECEF  Coordenadas geodesicas (WGS-84) -> ECEF [km]. Vectorizable.
%   lat_deg, lon_deg : latitud y longitud [deg] (pueden ser vectores)
%   h_km             : altura sobre el elipsoide [km] (escalar o vector)
%   r                : [M x 3] posiciones ECEF

a  = cfg.const.Re;
e2 = cfg.const.e2;
lat = deg2rad(lat_deg(:));
lon = deg2rad(lon_deg(:));

Nn = a ./ sqrt(1 - e2*sin(lat).^2);          % radio de curvatura primer vertical
X = (Nn + h_km).*cos(lat).*cos(lon);
Y = (Nn + h_km).*cos(lat).*sin(lon);
Z = (Nn.*(1 - e2) + h_km).*sin(lat);

r = [X, Y, Z];
end
