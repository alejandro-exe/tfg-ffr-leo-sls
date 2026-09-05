%% RUN_INTERFERENCE_DEMO  Interferencia co-canal (reuso-1): densidad x apuntamiento.
%  Anticipo de E3: ¿como degrada la SINR al saturar la constelacion, con patrones
%  de antena REALISTAS (terminal ITU-R S.1428 + haz satelite Bessel con suelo)?
%  Barre densidad de satelites y modo de apuntamiento del interferente.
%
%  METRICA de H1 (clave): la PENALIZACION POR INTERFERENCIA = C/N - SINR (el
%  "hueco" reuso-1). No es la SINR absoluta (que puede subir con la densidad por
%  mejor geometria del servidor), sino cuanto le quita la interferencia.
%
%  Modos de apuntamiento del interferente:
%    'nadir'     -> haz a su sub-satellite point (NOMINAL).
%    'boresight' -> haz al usuario, Gmax (COTA SUPERIOR, peor caso).
%
%  Validacion esperada:
%    - SINR <= C/N siempre; penalizacion >= 0.
%    - la penalizacion CRECE con la densidad (mas interferentes co-canal).
%    - penalizacion(boresight) >= penalizacion(nadir) (cota superior).
%
%  Uso: ejecutar en la raiz (reconstruye geometria para cada caso).

clear; clc; close all;

FIGDIR = 'figs_interf';      % PNG a 300 dpi via save_fig (fuente unica)

%% Barrido de densidad (Walker i:T/P) y modos de apuntamiento
Tlist = [66, 300, 600, 1584];
Plist = [ 6,  12,  24,   72];
modes = {'nadir','boresight'};
nD = numel(Tlist);  nM = numel(modes);

sinr_mean = nan(nM,nD);  sinr_p5 = nan(nM,nD);  cn_mean = nan(nM,nD);
pen_mean  = nan(nM,nD);  pen_p5  = nan(nM,nD);  nint_mean = nan(nM,nD);
store = cell(nM,nD);

fprintf('===== BARRIDO DENSIDAD x APUNTAMIENTO (interferencia reuso-1, usuario 1) =====\n');
fprintf('Patrones: terminal RX = ITU-R S.1428 ; haz satelite = Bessel + suelo %d dB\n\n', ...
    config_default().radio.sidelobe_floor_dB);

for im = 1:nM
    fprintf('--- Apuntamiento interferente: %s ---\n', modes{im});
    for d = 1:nD
        cfg = config_default();
        cfg.constellations(1).T = Tlist(d);
        cfg.constellations(1).P = Plist(d);
        cfg.interf.satPointing  = modes{im};

        tvec   = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
        sats   = build_constellation(cfg);
        R_eci  = propagate(cfg, sats, tvec);
        R_ecef = eci2ecef(cfg, R_eci, tvec);
        users  = build_user_grid(cfg);
        G      = compute_geometry(cfg, users, R_ecef);
        S      = associate_serving(cfg, G);
        L      = compute_link_budget(cfg, S);
        INTF   = compute_interference(cfg, sats, users, G, S, L);

        sinr = INTF.SINR_dB(1,:);  cn = INTF.CN_dB(1,:);
        ni   = INTF.nInterf(1,:);  pen = cn - sinr;            % penalizacion >= 0
        ok   = ~isnan(sinr);

        sinr_mean(im,d) = mean(sinr(ok));   sinr_p5(im,d) = prctile(sinr(ok),5);
        cn_mean(im,d)   = mean(cn(ok));
        pen_mean(im,d)  = mean(pen(ok));    pen_p5(im,d)  = prctile(pen(ok),5);
        nint_mean(im,d) = mean(ni(ok));
        store{im,d} = struct('tvec',tvec,'sinr',sinr,'cn',cn,'pen',pen,'ok',ok);

        fprintf(['T=%4d (P=%2d): SINR med=%6.2f p5=%6.2f | C/N med=%6.2f | ' ...
                 'PENAL med=%5.2f p5=%5.2f dB | nInt=%.2f\n'], ...
            Tlist(d), Plist(d), sinr_mean(im,d), sinr_p5(im,d), cn_mean(im,d), ...
            pen_mean(im,d), pen_p5(im,d), nint_mean(im,d));
    end
    fprintf('\n');
end

fprintf('================================================================\n');
for im = 1:nM
    fprintf('PENALIZACION p5 @ T=1584 (%s): %.2f dB   (media %.2f dB)\n', ...
        modes{im}, pen_p5(im,end), pen_mean(im,end));
end
fprintf('Interferentes medios @ T=1584: %.2f\n', nint_mean(1,end));
fprintf('================================================================\n\n');

%% Figura (a) - Penalizacion por interferencia vs densidad (debe CRECER)
f1 = figure('Name','Penalizacion vs densidad','Color','w');
plot(Tlist, pen_mean(1,:), 'o-',  'LineWidth',1.6, 'DisplayName','nadir (nominal) - media'); hold on;
plot(Tlist, pen_mean(2,:), 's-',  'LineWidth',1.6, 'DisplayName','boresight (cota) - media');
plot(Tlist, pen_p5(2,:),   's--', 'LineWidth',1.2, 'DisplayName','boresight - p5');
grid on; xlabel('N.o de satelites en la constelacion'); ylabel('Penalizacion C/N - SINR [dB]');
title('Penalizacion por interferencia co-canal vs densidad (reuso-1)');
legend('Location','northwest');
save_fig(f1, FIGDIR, 'interf_a_penalizacion_vs_densidad');

%% Figura (b) - SINR p5 vs densidad, banda nominal..cota
f2 = figure('Name','SINR p5: banda nominal-cota','Color','w');
xx = [Tlist, fliplr(Tlist)];
yy = [sinr_p5(2,:), fliplr(sinr_p5(1,:))];     % entre boresight (peor) y nadir (mejor)
fill(xx, yy, [.85 .9 1], 'EdgeColor','none', 'FaceAlpha',0.6, ...
    'DisplayName','Banda nominal..cota'); hold on;
plot(Tlist, sinr_p5(1,:), 'o-', 'LineWidth',1.6, 'DisplayName','SINR p5 - nadir (nominal)');
plot(Tlist, sinr_p5(2,:), 's-', 'LineWidth',1.6, 'DisplayName','SINR p5 - boresight (cota)');
grid on; xlabel('N.o de satelites en la constelacion'); ylabel('SINR percentil 5 [dB]');
title('SINR de borde (p5) vs densidad: nominal vs cota superior');
legend('Location','best');
save_fig(f2, FIGDIR, 'interf_b_sinr_p5_banda_nadir_boresight');

%% Figura (c) - SINR vs C/N en el tiempo a densidad maxima, modo boresight
sb = store{2,nD};  tv = sb.tvec/60;
f3 = figure('Name','SINR vs C/N (boresight, T max)','Color','w');
plot(tv, sb.cn,   'b-', 'LineWidth',1.4, 'DisplayName','C/N (sin interferencia)'); hold on;
plot(tv, sb.sinr, 'r-', 'LineWidth',1.4, 'DisplayName','SINR (boresight)');
tt = tv(sb.ok); cc = sb.cn(sb.ok); ss = sb.sinr(sb.ok);
fill([tt fliplr(tt)], [cc fliplr(ss)], [1 .8 .8], 'EdgeColor','none', ...
    'FaceAlpha',0.5, 'DisplayName','Hueco de interferencia');
grid on; xlabel('Tiempo [min]'); ylabel('[dB]');
title(sprintf('SINR vs C/N (T=%d, boresight, usuario 1)', Tlist(nD)));
legend('Location','southeast');
save_fig(f3, FIGDIR, 'interf_c_hueco_sinr_vs_cn_boresight');

fprintf('\n[figuras] 3 PNG a 300 dpi en %s%s\n', FIGDIR, filesep);

%% Guardar resultados
save('interference_results.mat','Tlist','Plist','modes', ...
     'sinr_mean','sinr_p5','cn_mean','pen_mean','pen_p5','nint_mean','store');
fprintf('Resultados guardados en interference_results.mat\n');
