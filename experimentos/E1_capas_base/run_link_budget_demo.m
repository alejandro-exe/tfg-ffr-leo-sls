%% RUN_LINK_BUDGET_DEMO  Demostracion y validacion del modulo de radioenlace.
%  Carga la geometria (geometry_results.mat), calcula el presupuesto de
%  enlace downlink Ku y muestra graficas de validacion.
%
%  Validacion esperada (ver apartado de resultados al final):
%    - C/N maximo cuando la elevacion es maxima (rango minimo)
%    - C/N cae cuando el satelite se acerca a la mascara de elevacion
%    - C/N lluvia < C/N claro (varios dB de diferencia en Ku)
%    - El margen puede ser negativo en pases rasantes con lluvia intensa
%
%  Uso: ejecutar tras run_geometry_demo (requiere geometry_results.mat)

clear; clc; close all;

FIGDIR = 'figs_link';        % PNG a 300 dpi via save_fig (fuente unica)

%% 0. Cargar geometria
if ~exist('geometry_results.mat','file')
    fprintf('geometry_results.mat no encontrado. Ejecutando geometria primero...\n');
    run_geometry_demo;
end
load('geometry_results.mat','cfg','tvec','S');

%% 1. Calcular presupuesto de enlace
fprintf('Calculando presupuesto de enlace...\n');
L = compute_link_budget(cfg, S);

if L.p618_used
    fprintf('  Modelo atmosferico: ITU-R P.618 via toolbox\n');
else
    fprintf('  Modelo atmosferico: aproximacion simplificada (toolbox no disponible)\n');
end

%% 2. Estadisticas de validacion
fprintf('\n===== BALANCE DE ENLACE (usuario 1) =====\n');
fprintf('Parametros:\n');
fprintf('  Frecuencia: %.0f GHz   |  Ancho de banda: %.0f MHz\n', ...
    cfg.radio.freq_GHz, cfg.radio.B_MHz);
fprintf('  EIRPdensidad: %.0f dBW/MHz  |  G/T terminal: %.0f dB/K\n', ...
    cfg.radio.EIRPdensity_dBWMHz, cfg.radio.GT_dBK);
fprintf('  C/N requerido: %.0f dB   |  Disponibilidad lluvia: %.1f%%\n\n', ...
    cfg.radio.CNreq_dB, 100 - cfg.radio.p618_availability);

% Extraer series temporales del usuario 1 (sin NaN)
cn_c  = L.CN_clear_dB(1,:);
cn_r  = L.CN_rain_dB(1,:);
mg_c  = L.margin_clear_dB(1,:);
mg_r  = L.margin_rain_dB(1,:);
el_s  = S.el(1,:);
fspl  = L.FSPL_dB(1,:);
latm_c = L.Latm_clear_dB(1,:);
latm_r = L.Latm_rain_dB(1,:);
cap   = L.Cshannon_Mbps(1,:);

ok = ~isnan(cn_c);
fprintf('CIELO CLARO (excedencia 50%%):\n');
fprintf('  C/N maximo  = %6.1f dB  (elevacion %.1f deg, rango %.0f km)\n', ...
    max(cn_c(ok)), el_s(cn_c == max(cn_c(ok))), S.range(1, cn_c == max(cn_c(ok))));
fprintf('  C/N minimo  = %6.1f dB  (elevacion %.1f deg, rango %.0f km)\n', ...
    min(cn_c(ok)), el_s(cn_c == min(cn_c(ok))), S.range(1, cn_c == min(cn_c(ok))));
fprintf('  Margen min  = %6.1f dB  (%s)\n', min(mg_c(ok)), signo(min(mg_c(ok))));
fprintf('  Capacidad Shannon max = %.0f Mbps\n\n', max(cap(ok)));

fprintf('CON LLUVIA (excedencia %.1f%% = %.1f%% disponibilidad):\n', ...
    cfg.radio.p618_availability, 100 - cfg.radio.p618_availability);
fprintf('  C/N maximo  = %6.1f dB\n', max(cn_r(ok)));
fprintf('  C/N minimo  = %6.1f dB\n', min(cn_r(ok)));
fprintf('  Margen min  = %6.1f dB  (%s)\n', min(mg_r(ok)), signo(min(mg_r(ok))));
fprintf('  Atenuacion lluvia max = %.1f dB  (elev minima)\n', ...
    max(latm_r(ok) - latm_c(ok)));

degradacion = min(cn_c(ok)) - min(cn_r(ok));
fprintf('\nDegradacion maxima claro->lluvia: %.1f dB\n', degradacion);
if min(mg_r(ok)) < 0
    n_neg = sum(mg_r(ok) < 0);
    fprintf('AVISO: margen negativo en %d instantes (%.1f%%) con lluvia.\n', ...
        n_neg, 100*n_neg/sum(ok));
else
    fprintf('Margen positivo en todos los instantes incluso con lluvia.\n');
end
fprintf('==========================================\n\n');

%% 3. Figura 1 - C/N vs tiempo (cielo claro vs lluvia)
f1 = figure('Name','C/N vs tiempo','Color','w');
plot(tvec/60, cn_c, 'b-',  'LineWidth', 1.4, 'DisplayName','Cielo claro (50%)'); hold on;
plot(tvec/60, cn_r, 'r--', 'LineWidth', 1.4, 'DisplayName', ...
    sprintf('Lluvia (%.1f%% excedencia)', cfg.radio.p618_availability));
yline(cfg.radio.CNreq_dB, 'k:', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('C/N req = %.0f dB', cfg.radio.CNreq_dB));
grid on; xlabel('Tiempo [min]'); ylabel('C/N [dB]');
title('Evolucion temporal del C/N (downlink Ku, usuario 1)');
legend('Location','southeast'); ylim([0 max(cn_c(ok))+5]);
save_fig(f1, FIGDIR, 'link_a_cn_vs_tiempo_claro_vs_lluvia');

%% 4. Figura 2 - Margen de enlace vs tiempo
f2 = figure('Name','Margen de enlace','Color','w');
plot(tvec/60, mg_c, 'b-',  'LineWidth', 1.4, 'DisplayName','Cielo claro'); hold on;
plot(tvec/60, mg_r, 'r--', 'LineWidth', 1.4, 'DisplayName','Lluvia');
yline(0, 'k-', 'LineWidth', 1.0, 'DisplayName','Limite viabilidad');
grid on; xlabel('Tiempo [min]'); ylabel('Margen [dB]');
title('Margen de enlace (C/N - C/N_{req})'); legend('Location','southeast');
save_fig(f2, FIGDIR, 'link_b_margen_de_guarda_vs_tiempo');

%% 5. Figura 3 - Perdidas atmosfericas P.618 vs elevacion (tabla precalculada)
el_tbl = (cfg.geom.minElev : 2 : 89);
% Reconstruir la tabla (la funcion interna no la expone; extraer del scatter)
el_sorted = sort(unique(round(el_s(ok))));  % elevaciones reales redondeadas
f3 = figure('Name','Perdidas atmosfericas vs elevacion','Color','w');
scatter(el_s(ok), latm_c(ok), 20, 'b', 'filled', 'DisplayName','Cielo claro'); hold on;
scatter(el_s(ok), latm_r(ok), 20, 'r', 'filled', 'DisplayName', ...
    sprintf('Lluvia %.1f%%', cfg.radio.p618_availability));
grid on; xlabel('Elevacion [deg]'); ylabel('Atenuacion [dB]');
title('Perdidas atmosfericas ITU-R P.618 vs elevacion (lat=40.3N, Ku)');
legend('Location','northeast');
xlim([cfg.geom.minElev 90]);
save_fig(f3, FIGDIR, 'link_c_perdidas_p618_vs_elevacion');

%% 6. Figura 4 - FSPL vs rango + contribuciones al balance
f4 = figure('Name','Contribuciones al balance','Color','w');
subplot(2,1,1);
scatter(S.range(1,ok), fspl(ok), 15, el_s(ok), 'filled');
colorbar; xlabel('Slant range [km]'); ylabel('FSPL [dB]');
title('FSPL vs rango (color = elevacion)'); grid on;

subplot(2,1,2);
rng_ok = S.range(1,ok);
lc_ok  = latm_c(ok);
lr_ok  = latm_r(ok);
scatter(rng_ok, lc_ok, 15, 'b', 'filled', 'DisplayName','Latm claro'); hold on;
scatter(rng_ok, lr_ok, 15, 'r', 'filled', 'DisplayName','Latm lluvia');
xlabel('Slant range [km]'); ylabel('Perdidas atmosfericas [dB]');
title('Atenuacion P.618 vs rango'); legend; grid on;
save_fig(f4, FIGDIR, 'link_d_fspl_y_atmosfera_vs_rango');

fprintf('\n[figuras] 4 PNG a 300 dpi en %s%s\n', FIGDIR, filesep);

%% Guardar resultado
save('link_budget_results.mat','cfg','tvec','S','L');
fprintf('Resultados guardados en link_budget_results.mat\n');

%% Subfuncion auxiliar local (MATLAB R2016b+: funciones locales al final del script)
function s = signo(x)
if x >= 0, s = 'enlace viable'; else, s = 'ENLACE EN FALLO'; end
end
