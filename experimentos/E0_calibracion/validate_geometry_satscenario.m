function V = validate_geometry_satscenario(cfg, jsat)
%VALIDATE_GEOMETRY_SATSCENARIO  Validacion independiente de la geometria propia.
%   V = validate_geometry_satscenario()          % usa config_calib, satelite auto
%   V = validate_geometry_satscenario(cfg, jsat)  % cfg y satelite concretos
%
%   Reconstruye la orbita de UN satelite con la Satellite Communications Toolbox
%   (satelliteScenario, propagador 'two-body-keplerian', los MISMOS elementos
%   orbitales que build_constellation) y compara la ELEVACION y el RANGO al punto
%   de tierra con los que produce nuestra geometria propia (propagate + eci2ecef
%   + compute_geometry). Validacion INDEPENDIENTE (toolbox != aportacion).
%
%   Alineacion de epoca: la toolbox usa la orientacion real de la Tierra (GMST)
%   en su StartTime, mientras nuestro modelo usa theta_g0 (=0 por defecto). Para
%   comparar en el mismo marco, se MIDE el GMST0 de la toolbox en t0 (rotacion
%   inertial->ECEF) y se fija cfg.prop.theta_g0 = GMST0 SOLO en esta validacion.
%   No se toca config_default ni config_calib.
%
%   Tolerancias objetivo: pocos km en rango, decimas de grado en elevacion.
%
%   Salida (struct V): max/rms de |dEl| [deg] y |dRange| [km], y series temporales.

if nargin < 1 || isempty(cfg),  cfg = config_calib();  end

% --- Comprobacion de disponibilidad de la toolbox ---
% satelliteScenario es una CLASE (classdef): usar which/exist(...,'class'),
% NO exist(...,'file') (que devuelve 2 solo para ficheros .m sueltos).
if isempty(which('satelliteScenario'))
    error('validate_geometry_satscenario:noToolbox', ...
        'satelliteScenario no disponible (Satellite Communications Toolbox).');
end

Re = cfg.const.Re;
sats = build_constellation(cfg);
tvec = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;

%% 1. Escenario de la toolbox y GMST0 (la epoca se alinea ANTES de elegir jsat)
%   ORDEN CRITICO: si jsat se eligiera con la theta_g0 por defecto y luego se
%   trazara con theta_g0 = GMST0, serian DOS CONVENCIONES DE EPOCA DISTINTAS.
%   Desplazar theta_g0 decenas de grados mueve la posicion relativa del punto de
%   tierra respecto a la orbita, asi que el satelite elegido como "servidor en su
%   maxima elevacion" puede quedar POR DEBAJO DEL HORIZONTE en toda la ventana una
%   vez alineado. Por eso se mide GMST0 primero y se selecciona ya con cfg2.
T0  = datetime(2020,1,1,0,0,0);
sc  = satelliteScenario(T0, T0 + seconds(cfg.time.duration), cfg.time.dt);
gs  = groundStation(sc, cfg.ground.point(1), cfg.ground.point(2), ...
    'Altitude', cfg.ground.alt*1e3);        % Altitude en metros

% GMST0 = rotacion inertial->ECEF de la toolbox en t0. Se mide con un satelite
% SONDA (cualquiera sirve) para poder elegir jsat ya en la epoca correcta.
satProbe   = add_kepler_sat(sc, sats, 1, 'probe');
GMST0_sel  = measure_gmst0(satProbe, T0);
fprintf('GMST0 de seleccion (sonda j=1): %.4f deg\n', GMST0_sel);

cfg2 = cfg;  cfg2.prop.theta_g0 = GMST0_sel;   % epoca alineada YA para elegir

%% 2. Satelite a validar: el servidor en su maxima elevacion, YA bajo cfg2
if nargin < 2 || isempty(jsat)
    R_eci0  = propagate(cfg2, sats, tvec);
    R_ecef0 = eci2ecef(cfg2, R_eci0, tvec);
    users0  = build_user_grid(cfg2);
    G0      = compute_geometry(cfg2, users0, R_ecef0);
    S0      = associate_serving(cfg2, G0);
    [~,ke]  = max(S0.el(1,:));
    jsat    = S0.idx(1,ke);
    if isnan(jsat), jsat = 1; end
end
fprintf('Validando satelite j=%d (a=%.1f km, inc=%.1f deg, raan=%.1f deg, M0=%.1f deg)\n', ...
    jsat, sats.a(jsat), rad2deg(sats.inc(jsat)), rad2deg(sats.raan(jsat)), rad2deg(sats.M0(jsat)));

% El satelite jsat en la toolbox (con los MISMOS elementos que build_constellation)
if jsat == 1
    sat = satProbe;                         % la sonda ya es el satelite pedido
else
    sat = add_kepler_sat(sc, sats, jsat, 'val');
end

% GMST0 se RE-MIDE sobre jsat para el trazado, que es la convencion original.
% No es redundante: la rotacion GCRF->ECEF de la toolbox NO es una rotacion-z
% pura (precesion/nutacion), luego este estimador de un solo angulo depende
% ligeramente de la DIRECCION del satelite con que se mide. Esa dispersion ES
% el tilt de marco (~0.06 deg) que explica los residuos, y usar la sonda para
% trazar lo inyectaria como desalineacion de epoca -- amplificada por a/h (~11.6)
% cerca del cenit. La sonda sirve para ELEGIR (donde 0.06 deg es irrelevante);
% para COMPARAR se usa el propio jsat.
GMST0_deg = measure_gmst0(sat, T0);
fprintf('GMST0 medido de la toolbox en t0: %.4f deg (se usa como theta_g0)\n', GMST0_deg);
fprintf('  dispersion sonda-jsat = %.4f deg  (= tilt de marco GCRF vs ECI de rotacion-z)\n', ...
    GMST0_deg - GMST0_sel);
cfg2.prop.theta_g0 = GMST0_deg;             % epoca del TRAZADO/COMPARACION

%% 3. Elevacion y rango de la TOOLBOX (referencia independiente)
%   aer(gs,sat) sin tiempo devuelve la serie en los instantes de muestreo del
%   escenario (SampleTime = cfg.time.dt), que coinciden con tvec.
[~, el_tb, rng_tb] = aer(gs, sat);          % el en deg, range en m (vectores)
el_tb  = el_tb(:);
rng_tb = rng_tb(:) / 1e3;                   % km
Nt = numel(tvec);
if numel(el_tb) ~= Nt
    error('validate_geometry_satscenario:len', ...
        'La toolbox devolvio %d muestras, esperaba %d.', numel(el_tb), Nt);
end

%% 4. Nuestra geometria con la epoca alineada (theta_g0 = GMST0, ya fijada en 1)
R_eci  = propagate(cfg2, sats, tvec);
R_ecef = eci2ecef(cfg2, R_eci, tvec);
users  = build_user_grid(cfg2);
Rj     = R_ecef(jsat,:,:);                  % 1 x 3 x Nt (solo el satelite jsat)
Gj     = compute_geometry(cfg2, users, Rj);
el_own  = squeeze(Gj.el(1,1,:));
rng_own = squeeze(Gj.range(1,1,:));

%% 5. Comparacion
dEl  = el_own  - el_tb;
dRng = rng_own - rng_tb;
V.jsat        = jsat;
V.GMST0_deg   = GMST0_deg;
V.el_own = el_own;  V.el_tb = el_tb;
V.rng_own = rng_own; V.rng_tb = rng_tb;
V.maxAbsEl_deg  = max(abs(dEl));
V.rmsEl_deg     = sqrt(mean(dEl.^2));
V.maxAbsRng_km  = max(abs(dRng));
V.rmsRng_km     = sqrt(mean(dRng.^2));
% Rango de elevacion del paso: comprobacion de que el satelite elegido es
% REALMENTE VISIBLE en la ventana (max por encima de la mascara del perfil).
V.elMin_deg     = min(el_own);
V.elMax_deg     = max(el_own);
V.minElev_deg   = cfg.geom.minElev;
V.visible       = V.elMax_deg > cfg.geom.minElev;

fprintf('\n---------- VALIDACION GEOMETRIA vs satelliteScenario ----------\n');
fprintf('Paso      : elevacion propia en [%.2f, %.2f] deg (mascara %.1f deg) -> %s\n', ...
    V.elMin_deg, V.elMax_deg, V.minElev_deg, ternstr(V.visible,'VISIBLE','NO visible'));
fprintf('Elevacion : max|dEl| = %.4f deg | rms = %.4f deg\n', V.maxAbsEl_deg, V.rmsEl_deg);
fprintf('Rango     : max|dRng|= %.3f km  | rms = %.3f km\n',   V.maxAbsRng_km, V.rmsRng_km);
% Tolerancias: elevacion en decimas de grado; rango en pocos km. El residuo de
% rango se concentra CERCA DEL HORIZONTE, donde el rayo es rasante y un desajuste
% angular de ~0.03 deg (precesion/nutacion GCRF<->nuestro ECI idealizado de una
% sola rotacion-z, no modelada) se amplifica a varios km de rango. Es esperado.
tolEl = 0.1; tolRng = 10;
if V.maxAbsEl_deg <= tolEl && V.maxAbsRng_km <= tolRng
    fprintf('RESULTADO: COINCIDE (tol %.2f deg / %.0f km; residuo de rango = amplificacion rasante del tilt de marco ~0.03 deg).\n', tolEl, tolRng);
else
    fprintf('RESULTADO: revisar (fuera de tol %.2f deg / %.0f km).\n', tolEl, tolRng);
end

%% 6. Figura de contraste
figure('Name','Validacion geometria (propia vs satelliteScenario)','Color','w');
tmin = tvec/60;
subplot(2,1,1);
plot(tmin, el_own,'b-','LineWidth',1.4,'DisplayName','propia'); hold on;
plot(tmin, el_tb,'r--','LineWidth',1.2,'DisplayName','satelliteScenario');
grid on; ylabel('Elevacion [deg]'); legend('Location','best');
title(sprintf('Satelite j=%d: geometria propia vs toolbox', jsat));
subplot(2,1,2);
plot(tmin, rng_own,'b-','LineWidth',1.4,'DisplayName','propia'); hold on;
plot(tmin, rng_tb,'r--','LineWidth',1.2,'DisplayName','satelliteScenario');
grid on; ylabel('Rango [km]'); xlabel('Tiempo [min]'); legend('Location','best');
end

% ----------------------------------------------------------------------------
function sat = add_kepler_sat(sc, sats, j, name)
%ADD_KEPLER_SAT  Anade al escenario el satelite j con los MISMOS elementos
%   orbitales que usa build_constellation (orbita circular kepleriana).
a_m    = sats.a(j) * 1e3;                   % m
ecc    = 0;                                 % orbita circular
inc_d  = rad2deg(sats.inc(j));
raan_d = rad2deg(sats.raan(j));
argp_d = 0;                                 % circular: argp = 0
nu_d   = rad2deg(sats.M0(j));               % e=0 -> true anom = arg de latitud = M0
sat = satellite(sc, a_m, ecc, inc_d, raan_d, argp_d, nu_d, ...
    'OrbitPropagator','two-body-keplerian', 'Name', name);
end

% ----------------------------------------------------------------------------
function g = measure_gmst0(sat, T0)
%MEASURE_GMST0  Rotacion inertial->ECEF de la toolbox en t0, en deg [-180,180].
ri = states(sat, T0, 'CoordinateFrame','inertial');   % 3x1 m
re = states(sat, T0, 'CoordinateFrame','ecef');       % 3x1 m
g  = atan2d(ri(2),ri(1)) - atan2d(re(2),re(1));
g  = mod(g + 180, 360) - 180;
end

% ----------------------------------------------------------------------------
function s = ternstr(cond, sTrue, sFalse)
if cond, s = sTrue; else, s = sFalse; end
end
