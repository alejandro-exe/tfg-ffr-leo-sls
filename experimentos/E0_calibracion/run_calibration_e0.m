%% RUN_CALIBRATION_E0  Escenario de CALIBRACION E0 (3GPP TR 38.821, Ka Set-1 LEO-600).
%  Ejecuta el pipeline COMPLETO con config_calib y produce la CDF de CIR/SINR de
%  los usuarios del HAZ CENTRAL de un cluster de 19 haces, siguiendo la
%  metodologia de 38.821 Sec. 6.1.3 (usuarios distribuidos en el haz central,
%  reuso-1, interferencia INTRA + INTER satelite).
%
%  DOBLE VIA:  Ku (config_default) = estudio  |  Ka Set-1 (config_calib) = E0.
%
%  Uso: ejecutar en la raiz del proyecto.
%
%  Salida por consola: media, mediana y percentiles 5/50/95 de la CIR (metrica
%  reina de 38.821) y de la SINR total. Bloque marcado para pegar los valores de
%  referencia de 38.821 Sec. 6.1.3.

clear; clc; close all;

% CRONOMETRO DE EXTREMO A EXTREMO. El warmup de la tabla
% P.618 se mide APARTE y se descuenta: OJO, E0 usa config_calib (Ka, 20 GHz,
% minElev=30) y por tanto una CLAVE DE CACHE DISTINTA de la del perfil Ku, luego
% construye SU PROPIA tabla aunque en la misma sesion ya se hubiera construido la de Ku.
tRun = tic;

%% 1. Configuracion de calibracion
cfg = config_calib();
fprintf('===== CALIBRACION E0 : %s =====\n', cfg.scenario);
fprintf('Ka DL: f=%.0f GHz | EIRPdens=%.1f dBW/MHz | Gmax=%.1f dBi | HPBW=%.4f deg\n', ...
    cfg.radio.freq_GHz, cfg.radio.EIRPdensity_dBWMHz, cfg.radio.Gmax_dBi, cfg.radio.beamwidth3dB_deg);
fprintf('VSAT: Grx_max=%.1f dBi | Tsys=%.1f K | G/T=%.2f dB/K | minElev=%d deg\n', ...
    cfg.radio.Grx_max_dBi, cfg.radio.Tsys_K, cfg.radio.GT_dBK, cfg.geom.minElev);

% Warmup EXPLICITO de la tabla P.618 (coste fijo, fuera del cronometro del
% experimento). No cambia ningun numero: la cache esta verificada bit a bit
% (test_atm_cache) y E0 sigue siendo bit-identica.
tW = tic;  atm_loss_dB(cfg, 45);  tWarm = toc(tW);
fprintf('[warmup] tabla P.618 (Ka) construida en %.1f s (coste fijo, se descuenta)\n', tWarm);
tExp = tic;

%% 2. Pipeline geometrico + radioenlace + interferencia inter-satelite
tvec   = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
sats   = build_constellation(cfg);
R_eci  = propagate(cfg, sats, tvec);
R_ecef = eci2ecef(cfg, R_eci, tvec);
users  = build_user_grid(cfg);
G      = compute_geometry(cfg, users, R_ecef);
S      = associate_serving(cfg, G);
L      = compute_link_budget(cfg, S);
INTF   = compute_interference(cfg, sats, users, G, S, L);

%% 3. Layout multihaz de 19 haces (separacion normativa 38.821 Sec. 6.1.1)
BL = build_beam_layout(cfg);
fprintf('\nLayout: %d haces | separacion inter-haz s=%.4f deg (3GPP TR 38.821 Sec.6.1.1)\n', ...
    BL.nBeams, BL.spacing_deg);

%% 4. Instante de evaluacion del enlace: mejor geometria (max elevacion, usuario 1)
el1 = S.el(1,:);  cov = ~isnan(el1);
if ~any(cov), error('run_calibration_e0:noCoverage','Sin cobertura del usuario 1.'); end
[~,ke] = max(el1);                       % instante de max elevacion del servidor
fprintf('Instante de enlace: k=%d, elevacion servidor=%.2f deg, rango=%.1f km\n', ...
    ke, S.el(1,ke), S.range(1,ke));

% Portadora del servidor en el centro de celda (Grel=0), ruido termico e
% interferencia inter-satelite agregada en ese instante:
EIRP_const  = cfg.radio.EIRPdensity_dBWMHz + 10*log10(cfg.radio.B_MHz);   % dBW/haz
C_dBW       = EIRP_const + cfg.radio.Grx_max_dBi - L.FSPL_dB(1,ke) - L.Latm_clear_dB(1,ke);
N_dBW       = -228.6 + 10*log10(cfg.radio.Tsys_K) + 10*log10(radio_band_Hz(cfg));
I_inter_dBW = INTF.I_dBW(1,ke);          % -Inf si no hay interferentes inter-sat
fprintf('C=%.2f dBW | N=%.2f dBW | C/N=%.2f dB | I_inter=%.2f dBW | nInter=%d\n', ...
    C_dBW, N_dBW, C_dBW-N_dBW, I_inter_dBW, INTF.nInterf(1,ke));

%% 5. CIR/SINR intra-satelite del haz central (usuarios distribuidos)
CIR = compute_beam_cir(cfg, BL, C_dBW, N_dBW, I_inter_dBW);

%% 6. Estadisticos de calibracion
st = CIR.stats;
fprintf('\n---------- CIR/SINR del HAZ CENTRAL (%d usuarios) ----------\n', cfg.beams.nUsers);
fprintf('CIR intra-satelite [dB]  : media=%.2f  mediana=%.2f  p5=%.2f  p50=%.2f  p95=%.2f\n', ...
    st.CIR_mean, st.CIR_median, st.CIR_p5, st.CIR_p50, st.CIR_p95);
if isfield(st,'SINR_mean')
    fprintf('SINR total (intra+inter+N): media=%.2f  mediana=%.2f  p5=%.2f  p50=%.2f  p95=%.2f\n', ...
        st.SINR_mean, st.SINR_median, st.SINR_p5, st.SINR_p50, st.SINR_p95);
end

%% 7. Comparacion con valores de referencia REALES de 3GPP TR 38.821
%  Referencia: TR 38.821 V16.0.0, Tabla 6.1.1.2-1 (Calibration results on DL
%  transmissions), fila Ka-band Set-1 LEO-600, reuso-1.
%
%  METRICA PRIMARIA: DL Geometry SIR  <->  nuestra CIR del haz central.
%    - Nuestra CIR es puro cociente de patrones de haz (portadora / interferencia
%      intra-satelite), es decir una SIR de geometria: comparacion SIR con SIR.
%  METRICA SECUNDARIA: DL Geometry SINR (comprobacion). Por la nota de la tabla:
%      Geometry SINR = -10*log10( I/C + N/C )   sobre el ancho de banda configurado.
%    - NO se fuerza la comparacion de SINR: nuestro nivel de RUIDO no esta calibrado
%      al detalle del 3GPP. Como en nuestro escenario el ruido es DESPRECIABLE
%      frente a la interferencia intra (N/C << I/C), nuestra SINR ~ SIR, luego la
%      comparacion valida es SIR vs SIR. La SINR se reporta solo como referencia.
ref.SIR_p5  = -3.0;   ref.SIR_p50  = -1.0;   ref.SIR_p95  = 1.2;   % DL Geometry SIR (PRIMARIA)
ref.SINR_p5 = -3.2;   ref.SINR_p50 = -1.2;   ref.SINR_p95 = 1.0;   % DL Geometry SINR (secundaria)
ref.source  = '3GPP TR 38.821 V16.0.0, Tabla 6.1.1.2-1, DL, Ka Set-1 LEO-600, reuso-1';

fprintf('\n========== CALIBRACION vs 3GPP TR 38.821 (Tabla 6.1.1.2-1) ==========\n');
fprintf('Fuente: %s\n', ref.source);

% --- PRIMARIA: nuestra CIR (SIR de geometria) vs DL Geometry SIR del 3GPP ---
fprintf('\n[PRIMARIA] DL Geometry SIR (3GPP) vs CIR del haz central (propia):\n');
fprintf('  3GPP SIR   : p5=%+.2f  p50=%+.2f  p95=%+.2f dB\n', ref.SIR_p5, ref.SIR_p50, ref.SIR_p95);
fprintf('  Propio CIR : p5=%+.2f  p50=%+.2f  p95=%+.2f dB\n', st.CIR_p5, st.CIR_p50, st.CIR_p95);
dSIR = [st.CIR_p5-ref.SIR_p5, st.CIR_p50-ref.SIR_p50, st.CIR_p95-ref.SIR_p95];
fprintf('  DELTA      : p5=%+.2f  p50=%+.2f  p95=%+.2f dB  (max|Delta|=%.2f)\n', ...
    dSIR(1), dSIR(2), dSIR(3), max(abs(dSIR)));
tolSIR = 1.0;   % dB (calibracion superada si max|Delta| <= ~1 dB)
if max(abs(dSIR)) <= tolSIR
    fprintf('  RESULTADO  : CALIBRACION SUPERADA (dentro de %.1f dB).\n', tolSIR);
else
    fprintf('  RESULTADO  : revisar (fuera de %.1f dB).\n', tolSIR);
end

% --- SECUNDARIA: DL Geometry SINR (solo referencia, ruido no calibrado) ---
if isfield(st,'SINR_p50')
    fprintf('\n[SECUNDARIA] DL Geometry SINR (3GPP) vs SINR total (propia, ruido no calibrado):\n');
    fprintf('  3GPP SINR  : p5=%+.2f  p50=%+.2f  p95=%+.2f dB\n', ref.SINR_p5, ref.SINR_p50, ref.SINR_p95);
    fprintf('  Propio SINR: p5=%+.2f  p50=%+.2f  p95=%+.2f dB\n', st.SINR_p5, st.SINR_p50, st.SINR_p95);
    fprintf('  DELTA      : p5=%+.2f  p50=%+.2f  p95=%+.2f dB\n', ...
        st.SINR_p5-ref.SINR_p5, st.SINR_p50-ref.SINR_p50, st.SINR_p95-ref.SINR_p95);
    fprintf('  NOTA: comparacion informativa; la valida es SIR vs SIR (N/C << I/C).\n');
end
fprintf('=====================================================================\n');

%% 8. Figuras: CDF de CIR y SINR + layout de haces
figure('Name','CDF CIR/SINR haz central (E0)','Color','w');
[fc,xc] = ecdf(CIR.CIR_dB);
plot(xc, fc, 'LineWidth',1.8, 'DisplayName','CIR intra-satelite'); hold on;
if ~isempty(CIR.SINR_dB)
    [fs,xs] = ecdf(CIR.SINR_dB);
    plot(xs, fs, 'LineWidth',1.8, 'DisplayName','SINR total');
end
grid on; xlabel('[dB]'); ylabel('CDF');
title(sprintf('CDF haz central - Calib E0 (Ka Set-1 LEO-600, %d haces)', BL.nBeams));
legend('Location','southeast');

figure('Name','Layout 19 haces + usuarios','Color','w');
th = linspace(0,2*pi,60);
for b = 1:BL.nBeams
    c = BL.offset_deg(b,:);
    R = cfg.radio.beamwidth3dB_deg/2;                 % circulo -3 dB
    plot(c(1)+R*cos(th), c(2)+R*sin(th), 'b-'); hold on;
    text(c(1), c(2), sprintf('%d',b), 'HorizontalAlignment','center', 'FontSize',8);
end
% usuarios del haz central (reconstruimos su nube con la misma semilla)
Rcell = cfg.beams.cellFrac*BL.spacing_deg; rng(0);
rr = Rcell*sqrt(rand(cfg.beams.nUsers,1)); tt = 2*pi*rand(cfg.beams.nUsers,1);
plot(rr.*cos(tt), rr.*sin(tt), '.', 'Color',[.85 .3 .3], 'MarkerSize',3);
axis equal; grid on; xlabel('offset x [deg]'); ylabel('offset y [deg]');
title('Cluster de 19 haces (circulos -3 dB) y usuarios del haz central');

%% 9. Guardar
timing = struct('warmup_p618_s', tWarm, 'experiment_s', toc(tExp), 'total_s', toc(tRun));
fprintf('\n[tiempos] experimento %.1f s (+ %.1f s de warmup P.618) | total %.1f s\n', ...
    timing.experiment_s, timing.warmup_p618_s, timing.total_s);

save('calibration_e0_results.mat','cfg','BL','CIR','C_dBW','N_dBW','I_inter_dBW','ke','timing');
fprintf('\nResultados guardados en calibration_e0_results.mat\n');
