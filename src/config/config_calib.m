function cfg = config_calib()
%CONFIG_CALIB  Perfil de CALIBRACION E0 (3GPP TR 38.821, Ka-band Set-1, LEO-600).
%   cfg = config_calib()
%
%   Parte de config_default() (banda Ku, escenario de ESTUDIO) y SOBRESCRIBE solo
%   lo necesario para reproducir el escenario de referencia del 3GPP:
%     - TR 38.821 V16.0.0, Tabla 6.1.1.1-1  (Ka-band, Satellite parameters Set-1, LEO-600)
%     - TR 38.821 V16.0.0, Tabla 6.1.1.1-3    (terminal VSAT, Ka)
%     - TR 38.821 V16.0.0, Tabla 6.1.3.2-1  (elevacion objetivo LEO)
%
%   DOBLE VIA (importante para la defensa):
%     Ku (config_default)  -> escenario de ESTUDIO (FFR, saturacion, KPIs propios).
%     Ka Set-1 (este perfil) -> CALIBRACION E0 contra valores publicados del 3GPP.
%   config_default.m NO se toca: este perfil solo redefine campos sobre una copia.
%
%   Uso: cfg = config_calib();  (luego run_calibration_e0).

cfg = config_default();                    % base Ku, se sobrescribe abajo

%% ----- Identificacion del escenario -----
cfg.scenario = 'calib_38821_KaSet1_LEO600';

%% ----- Constelacion: LEO-600 (Tabla 6.1.1.1-1) -----
cfg.constellations(1).h   = 600;           % km  altitud LEO-600 de la tabla
cfg.constellations(1).inc = 53;            % deg  38.821 no fija el plano: mantenemos 53

%% ----- Radioenlace: Ka DL (Tabla 6.1.1.1-1, Set-1) -----
cfg.radio.freq_GHz           = 20;         % GHz  portadora DL Ka
% Densidad EIRP del satelite (Nota 5: incluye back-off de 5 dB). El valor de la
% tabla para Ka Set-1 LEO-600 es 4 dBW/MHz.
cfg.radio.EIRPdensity_dBWMHz = 4;          % dBW/MHz  (Tabla 6.1.1.1-1, Ka Set-1 LEO-600)
cfg.radio.Gmax_dBi           = 38.5;       % dBi  Satellite Tx max Gain (Set-1)
cfg.radio.beamwidth3dB_deg   = 1.7647;     % deg  3 dB beamwidth de la tabla (Set-1)

%% ----- Terminal VSAT en tierra (Tabla 6.1.1.1-3, Ka) -----
cfg.radio.Grx_max_dBi = 39.7;              % dBi  VSAT Rx antenna gain
% G/T derivado de T_antena y figura de ruido (Tabla 6.1.1.1-3):
%   Tsys = Tant + T0*(F - 1),   F = 10^(NF/10)
%   G/T  = Grx_max - 10*log10(Tsys)
Tant  = 150;                               % K   temperatura de antena (cielo claro, VSAT)
NF_dB = 1.2;                               % dB  figura de ruido del receptor
T0    = 290;                               % K   temperatura de referencia
cfg.radio.Tsys_K = Tant + T0*(10^(NF_dB/10) - 1);                 % ~242.3 K
cfg.radio.GT_dBK = cfg.radio.Grx_max_dBi - 10*log10(cfg.radio.Tsys_K);  % ~15.9 dB/K

%% ----- Geometria: elevacion objetivo LEO (Tabla 6.1.3.2-1) -----
cfg.geom.minElev = 30;                     % deg  elevacion minima de servicio (LEO)

%% ----- Layout multihaz (Sec. 6.1.1: 19 haces) -----
% La separacion angular entre haces la calcula build_beam_layout con la formula
% normativa  s_deg = rad2deg(sqrt(3)*sin(HPBW/2))  (NO "1.0 anchos de haz").
cfg.beams.nBeams   = 19;                   % 19 haces (centro + 2 anillos hex)
cfg.beams.nUsers   = 1000;                 % usuarios distribuidos en el haz central (CDF)
cfg.beams.cellFrac = 0.5;                  % radio de celda = cellFrac * s_deg (crossover en s/2)

%% ----- Nota atmosfera -----
% Cielo claro con A_cenit/sin(elev) y margen de sombreado 0 (VSAT), ya cubierto
% por atm_loss_dB (P.618). No se toca aqui.

end
