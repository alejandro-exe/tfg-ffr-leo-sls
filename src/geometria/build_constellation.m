function sats = build_constellation(cfg)
%BUILD_CONSTELLATION  Genera los elementos orbitales (Walker) de todos los satelites.
%   Recorre cfg.constellations (1 o varias) y devuelve un struct con vectores
%   columna de longitud N (numero total de satelites):
%       sats.a     semieje mayor [km]        (orbita circular: a = Re + h)
%       sats.inc   inclinacion [rad]
%       sats.raan  ascension recta del nodo [rad]
%       sats.M0    anomalia media inicial [rad]
%       sats.n     movimiento medio [rad/s]
%       sats.cid   indice de constelacion (1..C) de cada satelite
%       sats.h     altitud [km]
%       sats.minElev  elevacion minima de OPERACION de cada satelite [deg]
%       sats.N     numero total de satelites
%
%   MASCARA POR CONSTELACION (sats.minElev). Cada operador declara su propia
%   elevacion minima de servicio y NO son iguales: Starlink opera desde 25 deg
%   mientras OneWeb no baja de 55 deg (nominal global; 45-50 en bajas latitudes)
%   [OneWeb Technical Narrative, p.13 Sec. A.4]. Con una unica mascara global se
%   dejarian entrar como interferentes satelites ajenos a elevaciones a las que NO
%   transmiten, sobreestimando la probabilidad de evento in-line de E4.
%   Se toma de cfg.constellations(c).minElev si existe; si no, del global
%   cfg.geom.minElev (retrocompatible: una cfg sin el campo se comporta como antes).
%
%   Walker  i:T/P/F  ->  S = T/P satelites por plano.
%     - RAAN del plano p:  Om0 + (span/P)*p ,  span = 360 (delta) o 180 (star)
%     - Anomalia media del sat s del plano p:
%           M0 + 2*pi*s/S  +  2*pi*F*p/T      (desfase entre planos = Walker F)

mu = cfg.const.mu;  Re = cfg.const.Re;
a_all=[]; inc_all=[]; raan_all=[]; M0_all=[]; n_all=[]; cid_all=[]; h_all=[];
mel_all=[];

for c = 1:numel(cfg.constellations)
    K = cfg.constellations(c);
    if mod(K.T, K.P) ~= 0
        error('build_constellation:Walker', ...
              'En %s: T (%d) debe ser multiplo de P (%d).', K.name, K.T, K.P);
    end
    S   = K.T / K.P;                 % satelites por plano
    a   = Re + K.h;                  % km
    n   = sqrt(mu / a^3);            % rad/s
    inc = deg2rad(K.inc);
    if strcmpi(K.pattern,'star'),  span = 180;  else,  span = 360;  end

    % Mascara de elevacion PROPIA de esta constelacion (ver cabecera)
    if isfield(K,'minElev') && ~isempty(K.minElev)
        mel = K.minElev;
    else
        mel = cfg.geom.minElev;
    end

    for p = 0:K.P-1
        raan = deg2rad(K.Om0) + deg2rad(span) * p / K.P;
        for s = 0:S-1
            M0 = deg2rad(K.M0) + 2*pi*s/S + 2*pi*K.F*p/K.T;
            a_all(end+1,1)    = a;          %#ok<AGROW>
            inc_all(end+1,1)  = inc;        %#ok<AGROW>
            raan_all(end+1,1) = raan;       %#ok<AGROW>
            M0_all(end+1,1)   = mod(M0,2*pi);%#ok<AGROW>
            n_all(end+1,1)    = n;          %#ok<AGROW>
            cid_all(end+1,1)  = c;          %#ok<AGROW>
            h_all(end+1,1)    = K.h;        %#ok<AGROW>
            mel_all(end+1,1)  = mel;        %#ok<AGROW>
        end
    end
end

sats.a    = a_all;   sats.inc = inc_all;  sats.raan = raan_all;
sats.M0   = M0_all;  sats.n   = n_all;    sats.cid  = cid_all;
sats.h    = h_all;   sats.N   = numel(a_all);
sats.minElev = mel_all;
end
