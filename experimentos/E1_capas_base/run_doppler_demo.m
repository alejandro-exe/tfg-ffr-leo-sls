%% RUN_DOPPLER_DEMO  Demostracion y validacion del modulo Doppler.
%  Carga la geometria (geometry_results.mat), recalcula posiciones/velocidades
%  ECI/ECEF (geometry_results NO guarda R_ecef ni V_eci) y obtiene el corrimiento
%  Doppler del satelite servidor sobre el usuario 1.
%
%  Validacion esperada:
%    - df cruza CERO en el instante de maxima elevacion (rdot=0, acercamiento min.)
%    - |df| maximo cerca de la mascara de elevacion (mayor velocidad radial)
%    - orden de magnitud: cientos de kHz en banda Ku (12 GHz, LEO 550 km)
%    - rdot analitico ~ diferencia finita de S.range (comprobacion cruzada)
%
%  Uso: ejecutar tras run_geometry_demo (requiere geometry_results.mat)

clear; clc; close all;

%% 0. Cargar geometria
if ~exist('geometry_results.mat','file')
    fprintf('geometry_results.mat no encontrado. Ejecutando geometria primero...\n');
    run_geometry_demo;
end
load('geometry_results.mat','cfg','tvec','sats','users','S');

%% 1. Recalcular posiciones y velocidades (geometry_results no las guarda)
[R_eci, V_eci] = propagate(cfg, sats, tvec);
R_ecef = eci2ecef(cfg, R_eci, tvec);

%% 2. Calcular Doppler del servidor
fprintf('Calculando Doppler del satelite servidor...\n');
dop = compute_doppler(cfg, sats, S, R_ecef, V_eci, tvec);

%% 3. Serie temporal del usuario 1
m  = 1;
df = dop.df_Hz(m,:);          % Hz
rd = dop.rdot_ms(m,:);        % m/s
el = S.el(m,:);               % deg
ok = ~isnan(df);

[~, kmax_el] = max(el);                          % instante de maxima elevacion
[absdf_max, k_absmax] = max(abs(df));            % instante de |df| maximo

fprintf('\n===== DOPPLER (usuario 1, servidor) =====\n');
fprintf('Frecuencia portadora: %.0f GHz\n', cfg.radio.freq_GHz);
fprintf('  |df| maximo  = %7.1f kHz  (elev %.1f deg)\n', ...
    absdf_max/1e3, el(k_absmax));
fprintf('  df @ max elevacion (%.1f deg) = %+.2f kHz   <- debe ~0\n', ...
    el(kmax_el), df(kmax_el)/1e3);
fprintf('  rango de df  = [%+.1f, %+.1f] kHz\n', min(df(ok))/1e3, max(df(ok))/1e3);
absrd_max = max(abs(rd(ok)));
fprintf('  |rdot| maximo = %.0f m/s  (%.2f km/s)\n', absrd_max, absrd_max/1e3);

%% 4. Comprobacion cruzada: rdot analitico vs diferencia finita de S.range
%    rdot ~ d(range)/dt. Solo es valido DENTRO de un mismo pase del servidor
%    (en los handover S.idx cambia y la diferencia finita salta).
rng_m = S.range(m,:) * 1e3;                       % km -> m
rdot_fd = nan(size(rng_m));
rdot_fd(2:end-1) = (rng_m(3:end) - rng_m(1:end-2)) ./ (tvec(3:end) - tvec(1:end-2));
no_ho = [false, diff(S.idx(m,:))==0];             % instantes sin handover previo
no_ho2 = no_ho & [no_ho(2:end), false];           % vecinos del mismo satelite
cmp = ok & no_ho2 & ~isnan(rdot_fd);
if any(cmp)
    err = rd(cmp) - rdot_fd(cmp);
    fprintf('  Validacion rdot (analitico vs dif. finita, mismo pase):\n');
    fprintf('     RMS error = %.2f m/s   max |error| = %.2f m/s\n', ...
        sqrt(mean(err.^2)), max(abs(err)));
end
fprintf('=========================================\n\n');

%% 5. Figura 1 - df Doppler vs tiempo
FIGDIR = 'figs_doppler';     % PNG a 300 dpi via save_fig (fuente unica)

f1 = figure('Name','Doppler vs tiempo','Color','w');
yyaxis left;
plot(tvec/60, df/1e3, 'b-', 'LineWidth', 1.4); hold on;
yline(0, 'k:', 'LineWidth', 1.0);
ylabel('Corrimiento Doppler \Deltaf [kHz]');
yyaxis right;
plot(tvec/60, el, 'r--', 'LineWidth', 1.0);
plot(tvec(kmax_el)/60, el(kmax_el), 'rv', 'MarkerFaceColor','r', 'MarkerSize',8);
ylabel('Elevacion del servidor [deg]'); ylim([0 90]);
grid on; xlabel('Tiempo [min]');
title('Doppler del servidor (\Deltaf=0 en maxima elevacion; |\Deltaf| max en la mascara)');
legend('\Deltaf [kHz]','\Deltaf=0','Elevacion','Max elevacion','Location','best');
save_fig(f1, FIGDIR, 'doppler_a_corrimiento_y_elevacion_vs_tiempo');

%% 6. Figura 2 - df vs elevacion (forma caracteristica del pase)
f2 = figure('Name','Doppler vs elevacion','Color','w');
scatter(el(ok), df(ok)/1e3, 18, tvec(ok)/60, 'filled');
cb = colorbar; cb.Label.String = 'Tiempo [min]';
grid on; xlabel('Elevacion del servidor [deg]'); ylabel('\Deltaf [kHz]');
title('Doppler vs elevacion (cruce por 0 en el cenit del pase)');
save_fig(f2, FIGDIR, 'doppler_b_corrimiento_vs_elevacion');

fprintf('\n[figuras] 2 PNG a 300 dpi en %s%s\n', FIGDIR, filesep);

%% Guardar resultado
save('doppler_results.mat','cfg','tvec','S','dop');
fprintf('Resultados guardados en doppler_results.mat\n');
