%% REGEN_CAP8   Figuras 8.1 a 8.4 + PANEL DOBLE 8.5 (Tarea 4), SIN TITULO.
%
%  Todas CATEGORIA A: lo graficado se LEE de minelev_results.mat / ffr_results.mat /
%  oneweb_results.mat. No se llama a compute_sinr_ffr, compute_interference,
%  compute_kpis, ffr_allocate, ffr_policy, run_one_density ni run_sweep_points.
%
%  TAREA 4 (figura 8.5): 'elevacion' y 'densificacion' cuentan como UNA figura de
%  panel doble. Se REGENERA desde oneweb_results.mat (no se compone a partir de los
%  dos PNG). Disposicion LADO A LADO: los dos paneles comparten la MISMA magnitud
%  en el eje Y (SINR_edge_p5) y solo cambian de eje X, asi que en horizontal se
%  comparan de un vistazo y la figura queda en proporcion ~2:1, apta para media
%  pagina; apilada quedaria alta y desperdiciaria el ancho de caja.
%
%  Uso: run('.../figuras/scripts/regen_cap8.m')

HERE     = fileparts(mfilename('fullpath'));
PROJROOT = fileparts(fileparts(HERE));
addpath(PROJROOT);
STAGE = fullfile(PROJROOT,'REDACCION','_staging');
FINAL = fullfile(PROJROOT,'REDACCION','figuras');
if ~exist(STAGE,'dir'), mkdir(STAGE); end
promote = @(n) copyfile(fullfile(STAGE,[n '.png']), fullfile(FINAL,[n '.png']));
bigfont = @(fh) set(findall(fh,'-property','FontSize'),'FontSize',11);
panel   = @(ax,t) text(ax, 0.02, 0.98, t, 'Units','normalized','FontSize',12, ...
            'FontWeight','bold','VerticalAlignment','top','HorizontalAlignment','left');
W2 = [100 100 760 360];      % 2 paneles
WL = [100 100 760 420];      % idem, con sitio para la leyenda compartida abajo

fprintf('\n################ CAPITULO 8 ################\n');

%% ===================== 8.1  Compresion angular (E3b) ====================
EB = load(fullfile(PROJROOT,'minelev_results.mat'), ...
          'elGrid','sepRad','sepTan','sepNom','STRAT','nm');
fprintf('\n--- 8.1 ---\n');
sn = interp1(EB.elGrid, EB.sepRad, [60 45 25 10]) / EB.sepNom;
fprintf('  sep/nom a 60/45/25/10 deg = %s  (publicado 0.771 0.528 0.210 0.053)\n', mat2str(round(sn,4)));
d81a = max(abs(sn - [0.771 0.528 0.210 0.053]));
fprintf('  max|dif| = %.5f (tol 0.0006)\n', d81a);
assert(d81a <= 0.0006, '8.1: la compresion analitica no reproduce la publicada.');
% Bins estratificados. Publicado: 10-20 -> 0.112 / -15.05 / -10.03 / -8.83
%                                 30-40 -> 0.356 /  -9.98 /  -4.28 / -2.47
%                                 60-90 -> 0.889 /  -4.35 /  +5.25 / +8.02
iB = [1 3 6];                               % bins centrados en 15, 35 y 75 deg
pubC = [0.112 0.356 0.889];
pubS = [-15.05 -9.98 -4.35;                 % reuse1   (col 1)
        -10.03 -4.28  5.25;                 % ffr(3)   (col 4)
         -8.83 -2.47  8.02];                % ffr(4)   (col 5)
gotS = [EB.STRAT.SINR(iB,1).'; EB.STRAT.SINR(iB,4).'; EB.STRAT.SINR(iB,5).'];
fprintf('  bins = %s | compr = %s\n', mat2str(EB.STRAT.bins), mat2str(round(EB.STRAT.compr,4)));
d81b = max(abs(EB.STRAT.compr(iB) - pubC));
d81c = max(max(abs(gotS - pubS)));
fprintf('  max|dif| compresion = %.5f | SINR estratificada = %.5f (tol 0.006)\n', d81b, d81c);
assert(d81b <= 0.006 && d81c <= 0.006, '8.1: los bins estratificados no coinciden.');
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off','Position',WL);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold on; grid on; box on;
plot(EB.elGrid, EB.sepRad, '-', 'LineWidth',1.8, 'DisplayName','radial (se escorza)');
plot(EB.elGrid, EB.sepTan, '--','LineWidth',1.4, 'DisplayName','transversal');
yline(EB.sepNom, ':','LineWidth',1.2,'HandleVisibility','off');
xlabel('elevacion del servidor [deg]'); ylabel('separacion efectiva [deg]');
legend('Location','southeast');
panel(ax1,'(a)');

ax2 = nexttile(tl);
ok = ~isnan(EB.STRAT.compr);
yyaxis left; hold on; grid on; box on;
plot(EB.STRAT.bins(ok), EB.STRAT.compr(ok), '-o','LineWidth',1.8,'MarkerSize',5, ...
     'DisplayName','compresion');
ylabel('separacion efectiva / nominal [-]');
yyaxis right;
% las dos curvas de SINR comparten el eje derecho y por tanto el COLOR: se
% distinguen ademas por estilo de linea y marcador, si no eran indistinguibles.
plot(EB.STRAT.bins(ok), EB.STRAT.SINR(ok,1), '-s','LineWidth',1.8,'MarkerSize',5, ...
     'DisplayName','SINR reuse1');
plot(EB.STRAT.bins(ok), EB.STRAT.SINR(ok,4), '--^','LineWidth',1.8,'MarkerSize',5, ...
     'DisplayName','SINR ffr(\Delta=3)');
ylabel('SINR_{edge,p5} [dB]');
xlabel('elevacion REAL del servidor [deg]');
panel(ax2,'(b)');

lgd = legend(ax2, 'Orientation','horizontal','NumColumns',3);
lgd.Layout.Tile = 'south';
bigfont(f); save_fig(f, STAGE, 'fig_8_01_compresion_angular'); close(f);
promote('fig_8_01_compresion_angular');

%% ===================== 8.2  Auto-similitud (OneWeb) =====================
%  [PENDIENTE-CAP8-AUTOSIMILITUD] La figura grafica la DESVIACION de G(vecino) frente
%  al asintotico, y ese "vecino" es el ADYACENTE (theta_off = s), primer co-canal SOLO
%  en reuso pleno. Con Delta=3 el primer co-canal esta a sqrt(3)*s (-17.954 dB) y con
%  Delta=4 a 2*s (-18.456 dB); la invariancia se cumple igual de bien en las tres
%  distancias (recorridos 0.010 / 0.014 / 0.028 dB), luego lo que hay que acotar al
%  redactar el capitulo es el NUMERO, no la propiedad. La etiqueta del eje se deja
%  como esta hasta fijar la formulacion definitiva.
OW = load(fullfile(PROJROOT,'oneweb_results.mat'), ...
          'AS','cfg','OUT','OUTR','nm','ELIST','RE','densB','RD');
fprintf('\n--- 8.2 ---\n');
fprintf('  HPBW = %s\n', mat2str(OW.AS.HPBW));
fprintf('  dev  = %s\n', mat2str(round(OW.AS.dev,4)));
pub82 = [0.000 0.001 0.005 0.027 0.107 0.253 0.417];
d82 = max(abs(OW.AS.dev - pub82));
fprintf('  max|dif| vs publicado = %.5f (tol 0.0006)\n', d82);
assert(d82 <= 0.0006, '8.2: la desviacion de auto-similitud no coincide.');
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off');
semilogx(OW.AS.HPBW, OW.AS.dev, 'o-','LineWidth',1.6); grid on; box on; hold on;
xline(2.08,'--','Starlink','FontSize',11,'LabelVerticalAlignment','bottom');
xline(OW.cfg.radio.beamwidth3dB_deg,'--','OneWeb','FontSize',11, ...
      'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
xlabel('ancho de haz a -3 dB [deg]');
ylabel('desviacion de G(vecino) vs asintotico [dB]');
bigfont(f); save_fig(f, STAGE, 'fig_8_02_autosimilitud'); close(f);
promote('fig_8_02_autosimilitud');

%% ===================== 8.3  Composicion del percentil 5 =================
E2 = load(fullfile(PROJROOT,'ffr_results.mat'), 'p5split','KK','cases');
nm2 = cellfun(@(c) c.name, E2.cases, 'UniformOutput', false);
compC = 100*[E2.p5split.fracC_R];
compB = 100 - compC;
base  = 100*E2.KK{1}.centerFrac;
fprintf('\n--- 8.3 ---\n');
fprintf('  %% de CENTRO bajo el p5 global = %s\n', mat2str(round(compC,2)));
fprintf('  linea base (centro en la poblacion) = %.2f %%\n', base);
% Publicado: "reuso-1 y reuso-D: lo fija el BORDE | FFR: lo fija el CENTRO (100%)"
assert(all(compC(1:3) < 10), '8.3: en reuse1/reuseD el cuello de botella deberia ser el BORDE.');
assert(all(abs(compC(4:5) - 100) < 1e-9), '8.3: en FFR el cuello de botella deberia ser el CENTRO (100%%).');
assert(abs(base - 50.09) <= 0.006, '8.3: centerFrac = %.4f, esperado 0.5009.', base/100);
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off');
hb = bar([compC(:) compB(:)], 'stacked');
hb(1).FaceColor = [0.85 0.33 0.10];
hb(2).FaceColor = [0.00 0.45 0.74];
hold on;
yline(base,'k--','LineWidth',1.6, ...
      'Label',sprintf('linea base: %.1f%% de centro en la poblacion', base), ...
      'FontSize',11,'LabelHorizontalAlignment','left');
set(gca,'XTickLabel',nm2,'XTickLabelRotation',20);
ylim([0 100]); grid on; box on;
ylabel('composicion bajo el percentil 5 global [%]');
legend({'usuarios de CENTRO','usuarios de BORDE'},'Location','east');
bigfont(f); save_fig(f, STAGE, 'fig_8_03_composicion_percentil5'); close(f);
promote('fig_8_03_composicion_percentil5');

%% ===================== 8.4  Topologia hexagonal vs 4x4 ==================
selT = 1:5;
vh = arrayfun(@(c) OW.OUT.K{c}.viab.SINR_edge_p5,  selT);
vr = arrayfun(@(c) OW.OUTR.K{c}.viab.SINR_edge_p5, selT);
fprintf('\n--- 8.4 ---\n');
fprintf('  hexagonal = %s\n', mat2str(round(vh,3)));
fprintf('  cuadrada  = %s\n', mat2str(round(vr,3)));
fprintf('  diferencia= %s\n', mat2str(round(vr-vh,3)));
% Publicado: hex reuse1 -3.98; D=3 +5.14 -> rect +1.90 (-3.23); D=4 +7.45 -> +9.07 (+1.62)
d84 = max(abs([vh(1) vh(2) vh(3) vr(2) vr(3)] - [-3.98 5.14 7.45 1.90 9.07]));
fprintf('  max|dif| vs publicado = %.5f (tol 0.006)\n', d84);
assert(d84 <= 0.006, '8.4: la tabla hex/rect no coincide con la publicada.');
assert(abs((vr(2)-vh(2)) - (-3.23)) <= 0.006 && abs((vr(3)-vh(3)) - 1.62) <= 0.006, ...
    '8.4: las diferencias -3.23 / +1.62 dB no se reproducen.');
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off');
hb = bar([vh(:) vr(:)]);
hb(1).FaceColor = [0.20 0.45 0.70];
hb(2).FaceColor = [0.85 0.45 0.15];
hold on;
yline(0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(-6.7,'r--','LineWidth',1.5,'HandleVisibility','off');
for i = 1:numel(selT)
    d = vr(i) - vh(i);
    text(i, max(vh(i),vr(i)) + 0.7, sprintf('%+.2f dB', d), ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',11, ...
        'Color', [0.6 0 0]*(d<0) + [0 0.45 0]*(d>=0));
end
set(gca,'XTickLabel',OW.nm(selT),'XTickLabelRotation',15);
grid on; box on; ylabel('SINR_{edge,p5} [dB]');
legend({'hexagonal (19 haces)','cuadrada 4x4 (16 haces)'},'Location','northwest');
bigfont(f); save_fig(f, STAGE, 'fig_8_04_topologia_hex_vs_rect'); close(f);
promote('fig_8_04_topologia_hex_vs_rect');

%% ============ 8.5  PANEL DOBLE (Tarea 4): elevacion + densificacion =====
selO = [1 2 4];                                   % reuse1, reuseD(3), ffr(3,a*)
sE = @(c) cellfun(@(R) R.K{c}.viab.SINR_edge_p5, OW.RE);
sD = @(c) cellfun(@(R) R.K{c}.viab.SINR_edge_p5, OW.RD);

fprintf('\n--- 8.5 (panel doble) ---\n');
fprintf('  ELIST = %s\n', mat2str(OW.ELIST));
fprintf('  (a) reuse1 vs mascara   = %s\n', mat2str(round(sE(1),3)));
fprintf('  densB = %s\n', mat2str(round(OW.densB,4)));
fprintf('  (b) reuse1 vs densidad  = %s\n', mat2str(round(sD(1),3)));
fprintf('  (b) ffr(3) vs densidad  = %s\n', mat2str(round(sD(4),3)));
% Publicado: recorrido del eje de elevacion (reuse1, 70->45) = 0.375 dB;
%            filas IDENTICAS por debajo de 55 deg (-3.994 a 50 y 45)
recE = max(sE(1)) - min(sE(1));
fprintf('  recorrido del eje de elevacion (reuse1) = %.4f dB (publicado 0.375)\n', recE);
assert(abs(recE - 0.375) <= 0.0006, '8.5a: el recorrido de elevacion no coincide.');
% densificacion publicada: reuse1 [-3.98 -4.07 -4.20 -4.30]; ffr(3) [+5.14 +5.05 +4.93 +4.62]
d85 = max([ max(abs(sD(1) - [-3.98 -4.07 -4.20 -4.30])), ...
            max(abs(sD(4) - [ 5.14  5.05  4.93  4.62])) ]);
fprintf('  max|dif| densificacion vs publicado = %.5f (tol 0.006)\n', d85);
assert(d85 <= 0.006, '8.5b: la tabla de densificacion no coincide.');
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off','Position',WL);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold on; grid on; box on;
for c = selO
    plot(OW.ELIST, sE(c), 'o-','LineWidth',1.6,'DisplayName',OW.nm{c});
end
yline(0,'k--','LineWidth',1.2,'HandleVisibility','off');
yline(-6.7,'r--','LineWidth',1.2,'HandleVisibility','off');
xline(55,':','nominal 55 deg','FontSize',11,'HandleVisibility','off', ...
      'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
set(gca,'XDir','reverse');
xlabel('mascara de elevacion [deg]'); ylabel('SINR_{edge,p5} [dB]');
panel(ax1,'(a)');

ax2 = nexttile(tl); hold on; grid on; box on;
for c = selO
    plot(OW.densB, sD(c), 'o-','LineWidth',1.6,'DisplayName',OW.nm{c});
end
yline(0,'k--','LineWidth',1.2,'HandleVisibility','off');
yline(-6.7,'r--','LineWidth',1.2,'HandleVisibility','off');
xlabel('densidad de haces [haces/1000 km^2]'); ylabel('SINR_{edge,p5} [dB]');
panel(ax2,'(b)');

lgd = legend(ax1, 'Orientation','horizontal','NumColumns',3);
lgd.Layout.Tile = 'south';
bigfont(f); save_fig(f, STAGE, 'fig_8_05_oneweb_panel'); close(f);
promote('fig_8_05_oneweb_panel');

fprintf('\n[cap 8] terminado.\n');
