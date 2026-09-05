function export_e0_figures()
%EXPORT_E0_FIGURES  Genera y guarda (PNG) todas las figuras del documento E0.
%   export_e0_figures()
%
%   Ejecuta el pipeline de calibracion (config_calib) y la validacion geometrica,
%   y exporta las figuras de la calibracion E0 a la carpeta figs_e0/.
%   Pensado para correr en modo -batch (figuras invisibles, exportgraphics).

outdir = fullfile(pwd,'figs_e0');
if ~exist(outdir,'dir'), mkdir(outdir); end
res = 300;                                   % dpi de exportacion (calidad de memoria;
                                             % mismo criterio que save_fig.m)

%% ---------- Pipeline de calibracion (identico a run_calibration_e0) ----------
cfg    = config_calib();
tvec   = cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration;
sats   = build_constellation(cfg);
R_eci  = propagate(cfg, sats, tvec);
R_ecef = eci2ecef(cfg, R_eci, tvec);
users  = build_user_grid(cfg);
G      = compute_geometry(cfg, users, R_ecef);
S      = associate_serving(cfg, G);
L      = compute_link_budget(cfg, S);
INTF   = compute_interference(cfg, sats, users, G, S, L);
BL     = build_beam_layout(cfg);

el1 = S.el(1,:);  [~,ke] = max(el1);
EIRP_const  = cfg.radio.EIRPdensity_dBWMHz + 10*log10(cfg.radio.B_MHz);
C_dBW       = EIRP_const + cfg.radio.Grx_max_dBi - L.FSPL_dB(1,ke) - L.Latm_clear_dB(1,ke);
N_dBW       = -228.6 + 10*log10(cfg.radio.Tsys_K) + 10*log10(radio_band_Hz(cfg));
I_inter_dBW = INTF.I_dBW(1,ke);
CIR         = compute_beam_cir(cfg, BL, C_dBW, N_dBW, I_inter_dBW);

% Valores de referencia 3GPP TR 38.821 V16.0.0, Tabla 6.1.1.2-1 (DL Geometry SIR,
% Ka Set-1 LEO-600, reuso-1) — mismos que en run_calibration_e0.m §7.
ref.SIR  = [-3.0, -1.0,  1.2];    % [p5 p50 p95] DL Geometry SIR (PRIMARIA)
ref.SINR = [-3.2, -1.2,  1.0];    % [p5 p50 p95] DL Geometry SINR (secundaria)
own.SIR  = [CIR.stats.CIR_p5,  CIR.stats.CIR_p50,  CIR.stats.CIR_p95];
own.SINR = [CIR.stats.SINR_p5, CIR.stats.SINR_p50, CIR.stats.SINR_p95];
pct = [5 50 95];

%% ---------- Fig 1: layout de 19 haces + usuarios del haz central ----------
f = figure('Visible','off','Color','w','Position',[100 100 620 560]);
th = linspace(0,2*pi,120);
Rb = cfg.radio.beamwidth3dB_deg/2;
for b = 1:BL.nBeams
    c = BL.offset_deg(b,:);
    col = [0.2 0.4 0.9]; if b==1, col=[0.85 0.2 0.2]; end
    plot(c(1)+Rb*cos(th), c(2)+Rb*sin(th), '-','Color',col,'LineWidth',1.3); hold on;
    text(c(1), c(2), sprintf('%d',b),'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end
Rcell = cfg.beams.cellFrac*BL.spacing_deg; rng(0);
rr = Rcell*sqrt(rand(cfg.beams.nUsers,1)); tt = 2*pi*rand(cfg.beams.nUsers,1);
plot(rr.*cos(tt), rr.*sin(tt), '.','Color',[0.9 0.5 0.2],'MarkerSize',3);
axis equal; grid on; box on;
xlabel('offset angular x [deg]'); ylabel('offset angular y [deg]');
title({'Cluster de 19 haces (circulos a -3 dB) y usuarios del haz central', ...
       sprintf('separacion s=%.3f deg (38.821 §6.1.1), HPBW=%.3f deg', BL.spacing_deg, cfg.radio.beamwidth3dB_deg)});
exportgraphics(f, fullfile(outdir,'fig1_beam_layout.png'),'Resolution',res); close(f);

%% ---------- Fig 2: CDF de CIR/SINR del haz central + referencia 3GPP ----------
f = figure('Visible','off','Color','w','Position',[100 100 680 480]);
[fc,xc] = ecdf(CIR.CIR_dB);  plot(xc,fc,'LineWidth',2,'DisplayName','CIR intra-satelite (propia)'); hold on;
[fs,xs] = ecdf(CIR.SINR_dB); plot(xs,fs,'LineWidth',2,'DisplayName','SINR total (propia)');
% Puntos de referencia 3GPP (DL Geometry SIR) en sus percentiles p5/p50/p95:
plot(ref.SIR, pct/100, 'ks', 'MarkerSize',9, 'MarkerFaceColor',[0.9 0.85 0.2], ...
     'LineWidth',1.2, 'DisplayName','3GPP TR 38.821 SIR (p5/p50/p95)');
for i=1:3
    text(ref.SIR(i)+0.06, pct(i)/100, sprintf('p%d',pct(i)), 'FontSize',8);
end
grid on; box on; xlabel('[dB]'); ylabel('CDF'); xlim([-3.2 2]);
title('CDF de CIR/SINR del haz central vs 3GPP TR 38.821 (Ka Set-1 LEO-600, reuso-1)');
legend('Location','southeast');
exportgraphics(f, fullfile(outdir,'fig2_cdf_cir_sinr.png'),'Resolution',res); close(f);

%% ---------- Fig 3: CIR vs offset radial del usuario (degradacion al borde) ----------
f = figure('Visible','off','Color','w','Position',[100 100 640 460]);
scatter(CIR.rho_deg, CIR.CIR_dB, 8, CIR.CIR_dB, 'filled'); hold on;
xline(Rcell,'k--','radio de celda (s/2)');
grid on; box on; colormap(parula); cb=colorbar; cb.Label.String='CIR [dB]';
xlabel('offset radial del usuario \rho [deg]'); ylabel('CIR intra-satelite [dB]');
title('CIR frente a la posicion del usuario en el haz central');
exportgraphics(f, fullfile(outdir,'fig3_cir_vs_rho.png'),'Resolution',res); close(f);

%% ---------- Fig 4: corte de patron (por que la CIR ~0 en reuso-1) ----------
f = figure('Visible','off','Color','w','Position',[100 100 640 460]);
s   = BL.spacing_deg; HPBW = cfg.radio.beamwidth3dB_deg; fl = cfg.radio.sidelobe_floor_dB;
xx  = linspace(-2*s, 2*s, 800);
gC  = beam_gain_dB(abs(xx),   HPBW, fl);      % haz central
gN  = beam_gain_dB(abs(xx-s), HPBW, fl);      % haz vecino (anillo 1)
plot(xx,gC,'-','LineWidth',1.8,'DisplayName','haz central (servidor)'); hold on;
plot(xx,gN,'-','LineWidth',1.8,'DisplayName','haz vecino (co-canal)');
gx = beam_gain_dB(s/2,HPBW,fl);
plot(s/2,gx,'ko','MarkerFaceColor','k','HandleVisibility','off');
text(s/2,gx+0.6,sprintf(' crossover %.1f dB',gx));
grid on; box on; ylim([max(fl,-25) 1]);
xlabel('offset angular [deg]'); ylabel('ganancia relativa [dB]');
title('Corte 1D del patron Bessel: haz servidor vs vecino co-canal');
legend('Location','south');
exportgraphics(f, fullfile(outdir,'fig4_beam_cut.png'),'Resolution',res); close(f);

%% ---------- Fig 5: geometria del servidor en el tiempo (contexto) ----------
f = figure('Visible','off','Color','w','Position',[100 100 640 460]);
tmin = tvec/60;
yyaxis left;  plot(tmin, S.el(1,:),'-','LineWidth',1.6); ylabel('Elevacion servidor [deg]');
yyaxis right; plot(tmin, S.range(1,:),'-','LineWidth',1.6); ylabel('Rango servidor [km]');
grid on; box on; xlabel('Tiempo [min]');
title(sprintf('Geometria del satelite servidor (usuario UC3M) — instante de enlace k=%d', ke));
exportgraphics(f, fullfile(outdir,'fig5_server_geometry.png'),'Resolution',res); close(f);

%% ---------- Fig 7: comparacion por percentil vs 3GPP (SIR) ----------
f = figure('Visible','off','Color','w','Position',[100 100 700 460]);
subplot(1,2,1);   % barras: 3GPP SIR vs propio CIR por percentil
b = bar([ref.SIR; own.SIR]', 'grouped');
b(1).FaceColor=[0.85 0.7 0.2]; b(2).FaceColor=[0.2 0.4 0.85];
set(gca,'XTickLabel',{'p5','p50','p95'}); grid on; box on;
ylabel('[dB]'); legend({'3GPP SIR','Propio CIR'},'Location','northwest');
title('DL Geometry SIR: 3GPP vs propio');
subplot(1,2,2);   % delta por percentil
dSIR = own.SIR - ref.SIR;
bd = bar(dSIR); bd.FaceColor=[0.3 0.6 0.3];
set(gca,'XTickLabel',{'p5','p50','p95'}); grid on; box on;
ylabel('\Delta = propio - 3GPP [dB]'); ylim([0 1.2]);
yline(1,'r--','tol ~1 dB');
for i=1:3, text(i, dSIR(i)+0.04, sprintf('%+.2f',dSIR(i)),'HorizontalAlignment','center','FontSize',9); end
title(sprintf('Desvio por percentil (max|\\Delta|=%.2f dB)', max(abs(dSIR))));
sgtitle('Calibracion vs 3GPP TR 38.821 Tabla 6.1.1.2-1 — SUPERADA','FontWeight','bold');
exportgraphics(f, fullfile(outdir,'fig7_calib_vs_3gpp.png'),'Resolution',res); close(f);

%% ---------- Fig 6: validacion geometrica vs satelliteScenario ----------
try
    V = validate_geometry_satscenario(cfg);
    close all;   % la validacion crea su propia figura visible; la descartamos
    f = figure('Visible','off','Color','w','Position',[100 100 660 620]);
    subplot(3,1,1);
    plot(tmin,V.el_own,'b-','LineWidth',1.5,'DisplayName','propia'); hold on;
    plot(tmin,V.el_tb,'r--','LineWidth',1.2,'DisplayName','satelliteScenario');
    grid on; box on; ylabel('Elevacion [deg]'); legend('Location','best');
    title(sprintf('Validacion geometrica — satelite j=%d', V.jsat));
    subplot(3,1,2);
    plot(tmin,V.rng_own,'b-','LineWidth',1.5,'DisplayName','propia'); hold on;
    plot(tmin,V.rng_tb,'r--','LineWidth',1.2,'DisplayName','satelliteScenario');
    grid on; box on; ylabel('Rango [km]'); legend('Location','best');
    subplot(3,1,3);
    plot(tmin,V.el_own-V.el_tb,'-','LineWidth',1.3,'DisplayName','\Delta elevacion [deg]'); hold on;
    plot(tmin,V.rng_own-V.rng_tb,'-','LineWidth',1.3,'DisplayName','\Delta rango [km]');
    grid on; box on; xlabel('Tiempo [min]'); ylabel('residuo'); legend('Location','best');
    exportgraphics(f, fullfile(outdir,'fig6_geom_validation.png'),'Resolution',res); close(f);
    fprintf('Validacion: max|dEl|=%.4f deg, max|dRng|=%.3f km\n', V.maxAbsEl_deg, V.maxAbsRng_km);
catch ME
    warning('export_e0_figures:val','No se pudo generar la figura de validacion: %s', ME.message);
end

fprintf('Figuras exportadas a %s\n', outdir);
fprintf('CIR: media=%.2f p5=%.2f p50=%.2f p95=%.2f dB | C/N=%.2f dB\n', ...
    CIR.stats.CIR_mean, CIR.stats.CIR_p5, CIR.stats.CIR_p50, CIR.stats.CIR_p95, C_dBW-N_dBW);
end
