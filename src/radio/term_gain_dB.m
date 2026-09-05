function Gdbi = term_gain_dB(cfg, phi_deg)
%TERM_GAIN_DB  Patron de ganancia del terminal RX en tierra (ITU-R S.1428-1).
%   Gdbi = term_gain_dB(cfg, phi_deg)
%
%   Envolvente del Recommendation ITU-R S.1428-1 para antenas de estacion
%   terrena FSS (reflector circular), en dBi ABSOLUTOS. phi_deg es el angulo
%   off-axis respecto al eje de punteria del terminal (apunta al servidor).
%   phi = 0 -> Grx_max (el enlace al servidor NO cambia; C/N intacto).
%
%   La toolbox NO provee este patron (comprobado: no existe s1428/earthStationAntenna),
%   por lo que se implementa la envolvente a mano. Aportacion: integracion del
%   patron de terminal en el balance de interferencia (la envolvente es estandar).
%
%   D/lambda se DERIVA de la ganancia pico y la eficiencia de apertura:
%       Gmax = 10*log10( eta * (pi * D/lambda)^2 )
%       D/lambda = sqrt( 10^(Gmax/10) / (eta * pi^2) )
%   Perfil Ku REAL (terminal de 60 cm a 12 GHz): Grx_max = 35.3 dBi, eta = 0.6
%   -> D/lambda ~ 23.9, coherente con el D/lambda geometrico 0.60/lambda = 24.0.
%   (El perfil anterior, Grx_max = 39 dBi, daba D/lambda ~ 36.6.)
%
%   Estructura de la envolvente (25 < D/lambda <= 100, ITU-R S.1428-1):
%       0      <= phi < phi_m  :  Gmax - 2.5e-3*(D/L * phi)^2   (lobulo principal)
%       phi_m  <= phi < phi_r  :  G1                            (meseta constante)
%       phi_r  <= phi < 33.1   :  29 - 25*log10(phi)            (envolvente)
%       33.1   <= phi < 80     :  -9  dBi                       (suelo)
%       80     <= phi <= 180   :  -4  dBi                       (suelo)
%   con  Gmax = Grx_max_dBi (forzado, phi=0 -> Grx_max),
%        G1   = 29 - 25*log10(95/(D/L)),
%        phi_m = (20/(D/L)) * sqrt(Gmax - G1),  phi_r = 95/(D/L).
%
%   REGIMEN 20 <= D/lambda <= 25 (lo dispara el perfil Ku REAL, D/L = 23.9):
%   S.1428-1 define para ese tramo la MISMA envolvente salvo en la COLA, donde da
%   -9 dBi para 120 < phi <= 180 en vez de -4 dBi. Este codigo mantiene la version
%   de 25<D/L<=100, luego para interferentes a mas de 120 deg del eje de punteria
%   devuelve 5 dB DE MAS: es CONSERVADOR (sobreestima la interferencia). El caso es
%   ALCANZABLE (dos satelites visibles con mascara de 25 deg pueden separarse hasta
%   180-2*25 = 130 deg), pero inmaterial mientras I_inter << N (~40 dB por debajo).
%   Se AVISA en vez de silenciarlo.
%
%   Acepta phi_deg escalar o matricial (vectorizado). NaN -> NaN.

Gmax = cfg.radio.Grx_max_dBi;
eta  = cfg.radio.term_efficiency;

% D/lambda derivado de la ganancia pico y la eficiencia
DL = sqrt( 10^(Gmax/10) / (eta * pi^2) );

% Aviso de REGIMEN, UNA SOLA VEZ por configuracion: esta funcion se llama dentro del
% bucle (usuario, instante) de compute_interference -- decenas de miles de veces por
% experimento -- y un warning por llamada inundaria la salida y ralentizaria el
% barrido. Se rearma solo si cambia (Gmax, eta), o a mano con `clear term_gain_dB`.
persistent warnedKey
key = sprintf('%.6g|%.6g', Gmax, eta);
if ~strcmp(warnedKey, key)
    warnedKey = key;
    if DL > 100
        warning('term_gain_dB:regime', ...
            ['D/lambda = %.1f > 100: envolvente S.1428 implementada para ' ...
             '20<=D/L<=100; resultado aproximado.'], DL);
    elseif DL <= 25
        warning('term_gain_dB:regime25', ...
            ['D/lambda = %.2f <= 25: se aplica la envolvente S.1428-1 de ' ...
             '25<D/L<=100. El tramo 20<=D/L<=25 solo difiere en phi>120 deg ' ...
             '(-9 en vez de -4 dBi), donde este codigo SOBREESTIMA la ' ...
             'interferencia en 5 dB (conservador).'], DL);
    end
end

% Parametros de la envolvente
G1    = 29 - 25*log10(95/DL);                 % nivel de la meseta
phi_m = (20/DL) * sqrt(max(Gmax - G1, 0));    % fin del lobulo principal
phi_r = 95/DL;                                % inicio de la envolvente 29-25log

phi = abs(phi_deg);                           % patron simetrico
Gdbi = nan(size(phi));

m1 =            phi <  phi_m;
m2 = phi >= phi_m & phi <  phi_r;
m3 = phi >= phi_r & phi <  33.1;
m4 = phi >= 33.1  & phi <  80;
m5 = phi >= 80    & phi <= 180;

Gdbi(m1) = Gmax - 2.5e-3 * (DL * phi(m1)).^2;
Gdbi(m2) = G1;
Gdbi(m3) = 29 - 25*log10(phi(m3));
Gdbi(m4) = -9;
Gdbi(m5) = -4;
end
