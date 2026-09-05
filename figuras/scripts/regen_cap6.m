%% REGEN_CAP6   Figuras del cap. 6 (6.1-6.7 y 6.9), SIN TITULO, a REDACCION/_staging
%  y promocion. Ademas emite por CONSOLA la tabla deducida de la Sec. 6.3.
%  La 6.8 antigua (fig_6_08_kpi_vs_operadores) queda como material de ANEXO: su PNG
%  se conserva en REDACCION/figuras pero YA NO se regenera aqui, porque su bloque lo
%  sustituye el de la 6.9 (ver esa seccion).
%
%  Reglas: sin title()/sgtitle(); paneles etiquetados (a)/(b); ejes con unidad;
%  ASCII; fuentes >= 11 pt; 300 dpi via save_fig.
%
%  REGLA DE RECALCULO: TODAS las figuras de este capitulo son de CATEGORIA A. Todo
%  lo graficado se LEE del .mat correspondiente; no se llama a compute_sinr_ffr,
%  compute_interference, compute_kpis, ffr_allocate, ffr_policy, run_one_density ni
%  run_sweep_points.
%
%  VERIFICACION NUMERICA: como los vectores graficados SON los almacenados, la
%  comparacion contra "lo que produjo el original" es una identidad (max|dif| = 0
%  por construccion). Lo que se comprueba de verdad aqui es algo INDEPENDIENTE: que
%  esos vectores reproducen las TABLAS PUBLICADAS del corpus. Cualquier desviacion
%  aborta y se reporta.
%
%  Uso: run('.../figuras/scripts/regen_cap6.m')

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
% En los mapas de veredictos las celdas llegan al borde y llevan texto dentro, asi
% que la etiqueta de panel va FUERA, justo encima del eje.
panelOut= @(ax,t) text(ax, 0.00, 1.03, t, 'Units','normalized','FontSize',12, ...
            'FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left');
W2 = [100 100 760 360];      % 2 paneles: 7.9 x 3.75 in
WL = [100 100 760 420];      % idem, con sitio para la leyenda compartida abajo
% Imprime el max|dif| medido y su tolerancia. El veredicto lo da el assert que
% sigue a cada llamada: si se supera la tolerancia, el script ABORTA y se reporta.
chk = @(tag,d,tol) fprintf('  %-36s max|dif| = %.5f   (tol %.4f)\n', tag, d, tol);

fprintf('\n################ CAPITULO 6 ################\n');

%% ============ 6.1 y 6.2  CDFs de E2 (ffr_results.mat) ====================
E2 = load(fullfile(PROJROOT,'ffr_results.mat'), 'KK','cases','cfg');
nC2 = numel(E2.cases);
nm2 = cellfun(@(c) c.name, E2.cases, 'UniformOutput', false);

% Publicado en la tabla de resultados de E2 (ventana completa, T=1584):
pub_Rp5_edge = [5.6 8.2 9.3 9.9 11.2];        % R_p5 BORDE [Mbps]
pub_SINRp5   = [-4.66 4.54 7.79 -3.57 -3.57]; % SINR_p5 GLOBAL [dB]
got_Rp5_edge = cellfun(@(K) K.edge.R_p5_Mbps, E2.KK);
got_SINRp5   = cellfun(@(K) K.all.SINR_p5,    E2.KK);

fprintf('\n--- 6.1 / 6.2 (E2) ---\n');
fprintf('  esquemas   : %s\n', strjoin(nm2,' | '));
fprintf('  R_p5 borde : %s  (publicado %s)\n', mat2str(round(got_Rp5_edge,2)), mat2str(pub_Rp5_edge));
fprintf('  SINR_p5    : %s  (publicado %s)\n', mat2str(round(got_SINRp5,2)),   mat2str(pub_SINRp5));
d61 = max(abs(got_Rp5_edge - pub_Rp5_edge));  chk('6.1 R_p5 borde vs publicado', d61, 0.05);
d62 = max(abs(got_SINRp5   - pub_SINRp5));    chk('6.2 SINR_p5 vs publicado',    d62, 0.005);
assert(d61 <= 0.05 && d62 <= 0.005, 'E2 no reproduce las tablas publicadas.');

% --- 6.1 CDF de throughput de BORDE ---
f = figure('Color','w','Visible','off');
for c = 1:nC2
    cdf = E2.KK{c}.edge.cdf_R;
    if all(isnan(cdf(:))), continue; end
    semilogx(cdf(:,1), cdf(:,2), 'LineWidth',1.6, 'DisplayName',nm2{c}); hold on;
end
grid on; box on; xlabel('throughput por usuario R [Mbps]'); ylabel('CDF [-]');
legend('Location','southeast');
bigfont(f); save_fig(f, STAGE, 'fig_6_01_cdf_throughput_borde'); close(f);
promote('fig_6_01_cdf_throughput_borde');

% --- 6.2 CDF de SINR ---
f = figure('Color','w','Visible','off');
for c = 1:nC2
    cdf = E2.KK{c}.all.cdf_SINR;
    plot(cdf(:,1), cdf(:,2), 'LineWidth',1.6, 'DisplayName',nm2{c}); hold on;
end
xline(E2.cfg.kpi.gamma0_dB,'k--','DisplayName',sprintf('\\gamma_0 = %g dB',E2.cfg.kpi.gamma0_dB));
grid on; box on; xlabel('SINR [dB]'); ylabel('CDF [-]');
legend('Location','southeast');
bigfont(f); save_fig(f, STAGE, 'fig_6_02_cdf_sinr_por_esquema'); close(f);
promote('fig_6_02_cdf_sinr_por_esquema');

%% ============ 6.3  Satelites visibles y KPI vs densidad orbital ==========
%  FUENTE NUEVA: e3_density_sweep.mat, el eje de densidad orbital DEL CAPITULO 6.
%  Antes la 6.3 salia de convergence_studyA.mat (FASE A), que barre el
%  MISMO eje T pero con otro proposito -- estudio de CONVERGENCIA Y MEMORIA: rejilla
%  de 6 km (M=137), Nt=61, dos esquemas y sin veredictos. La Sec. 6.2 de la memoria
%  ya no se apoya en ese barrido sino en e3_density_sweep (7 puntos, P escalado con
%  T, M=553, Nt=121, los cinco esquemas del capitulo mas la adaptativa), asi que la
%  figura tiene que salir del MISMO .mat que la tabla o el apartado queda con dos
%  fuentes. La version de FASE A no se pierde: se sigue regenerando mas abajo con
%  nombre de ARCHIVO, porque es la figura correcta del precedente de dimensionado.
DS = load(fullfile(PROJROOT,'e3_density_sweep.mat'), ...
          'results','Tlist','Plist','schemeCases','cfg');

fprintf('\n--- 6.3 (eje de densidad orbital, e3_density_sweep.mat) ---\n');
fprintf('  Tlist = %s\n', mat2str(DS.Tlist));
fprintf('  Plist = %s\n', mat2str(DS.Plist));

% GUARDA DE FUENTE (no es paranoia: existe una tanda archivada con P = 6 planos
% FIJOS que tiene los MISMOS nombres de variable y la MISMA Tlist, y solo se
% distingue por P).
% Cargarlo por error daria una figura perfectamente plausible pero con reuse1 a
% -11.07 dB en T=1584 en vez de -4.97, porque con 6 planos fijos el servidor queda
% sistematicamente bajo. Se comprueba P, que es lo unico que los separa.
assert(isequal(DS.Tlist, [66 300 600 1000 1584 2500 4000]), ...
    '6.3: Tlist inesperada (%s).', mat2str(DS.Tlist));
assert(isequal(DS.Plist, [3 12 24 40 72 100 160]), ...
    ['6.3: Plist = %s. Esto NO es el barrido nuevo (P escalado con T): parece el ' ...
     'archivado e3_density_sweep_P6FIJO.mat, que NO debe citarse.'], mat2str(DS.Plist));

nmD = cellfun(@(c) c.name, DS.schemeCases, 'UniformOutput', false);
nSd = numel(nmD);   nTd = numel(DS.Tlist);
need = {'reuse1','reuseD (D=3)','reuseD (D=4)', ...
        'FFR (D=3, a=0.4)','FFR (D=4, a=0.4)','FFR-ADAPTATIVA (D=3)'};
miss = need(~ismember(need, nmD));
assert(isempty(miss), '6.3: faltan esquemas en el .mat: %s', strjoin(miss,', '));
idxD = @(s) find(strcmp(nmD, s), 1);
iR1  = idxD('reuse1');
i3   = [idxD('reuseD (D=3)') idxD('FFR (D=3, a=0.4)') idxD('FFR-ADAPTATIVA (D=3)')];
i4   = [idxD('reuseD (D=4)') idxD('FFR (D=4, a=0.4)')];

% Todo lo graficado se LEE del .mat: ninguna cifra escrita a mano en el dibujo.
sinrD = nan(nSd, nTd);
for c = 1:nSd
    for i = 1:nTd, sinrD(c,i) = DS.results{i}.K{c}.viab.SINR_edge_p5; end   % K = ventana completa
end
nVis = cellfun(@(r) r.diag.nVis_mean, DS.results);
covF = cellfun(@(r) r.diag.covFrac,   DS.results);
gth6 = DS.cfg.viab.gamma_th_dB;   gfl6 = DS.cfg.viab.gamma_floor_dB;

fprintf('  nVis_mean = %s\n', mat2str(round(nVis,4)));
fprintf('  covFrac   = %s\n', mat2str(round(covF,4)));
fprintf('  reuse1    = %s\n', mat2str(round(sinrD(iR1,:),4)));
fprintf('  D=3       = %s\n', mat2str(round(sinrD(i3(1),:),4)));
fprintf('  D=4       = %s\n', mat2str(round(sinrD(i4(1),:),4)));

% (a) los cuatro vectores que se dibujan, contra el barrido publicado. Tolerancia =
%     redondeo a 4 decimales (5e-5) con un margen minimo.
ref_nVis = [ 0.8373   2.4420   4.8770   7.9177  12.5097  19.7087  31.5780];
ref_r1   = [-11.6127 -10.2120  -6.4775  -5.5516  -4.9694  -4.6104  -4.5088];
ref_d3   = [ -6.2185  -4.5633   0.8221   2.5669   3.7957   4.5584   4.7448];
ref_d4   = [ -4.6256  -2.7941   3.8207   6.0411   7.2539   7.7355   7.8594];
tol4 = 6e-5;
d63a = max(abs(nVis            - ref_nVis)); chk('6.3 nVis_mean vs barrido',  d63a, tol4);
d63b = max(abs(sinrD(iR1,:)    - ref_r1  )); chk('6.3 SINR_p5 reuse1',        d63b, tol4);
d63c = max(abs(sinrD(i3(1),:)  - ref_d3  )); chk('6.3 SINR_p5 D=3',           d63c, tol4);
d63d = max(abs(sinrD(i4(1),:)  - ref_d4  )); chk('6.3 SINR_p5 D=4',           d63d, tol4);
assert(max([d63a d63b d63c d63d]) <= tol4, ...
    '6.3 no reproduce el barrido de densidad orbital publicado.');

% (b) COLAPSO POR DELTA -- se MIDE antes de colapsar. Las tres curvas de Delta=3
%     (reuso-3, FFR(3) y FFR adaptativa) y las dos de Delta=4 coinciden porque en el
%     BORDE la FFR ve el mismo conjunto co-canal que reuso-Delta (hallazgo 7 del
%     TFG). Dibujar UNA curva por Delta es por tanto correcto, pero si alguna vez
%     dejaran de coincidir seria un RESULTADO y no un detalle de dibujo: por eso el
%     assert, en lugar de colapsar a ciegas.
colD3 = max(max(abs(sinrD(i3,:) - sinrD(i3(1),:))));
colD4 = max(max(abs(sinrD(i4,:) - sinrD(i4(1),:))));
chk('6.3 colapso D=3 (3 esquemas)', colD3, 5e-5);
chk('6.3 colapso D=4 (2 esquemas)', colD4, 5e-5);
assert(colD3 <= 5e-5 && colD4 <= 5e-5, ...
    ['6.3: los esquemas del mismo Delta YA NO coinciden (D=3: %.3e, D=4: %.3e dB). ' ...
     'No colapsar las curvas: es un resultado.'], colD3, colD4);

% (c) ANCLA. La fila T=1584 debe reproducir ffr_results.mat (E2), el escenario de
%     cabecera de la Sec. 6.1. run_density_sweep ya lo verifica al ejecutarse (o no
%     guarda el .mat); aqui se REPITE de forma independiente, leyendo los dos .mat y
%     emparejando por NOMBRE, sobre los tres valores que la figura dibuja.
E2a  = load(fullfile(PROJROOT,'ffr_results.mat'), 'cases','KK');
nm2a = cellfun(@(c) c.name, E2a.cases, 'UniformOutput', false);
iNom = find(DS.Tlist == 1584, 1);
dAnc = 0;
for pair = {{'reuse1',iR1}, {'reuseD (D=3)',i3(1)}, {'reuseD (D=4)',i4(1)}}
    j = find(strcmp(nm2a, pair{1}{1}), 1);
    assert(~isempty(j), '6.3: "%s" no esta en ffr_results.mat.', pair{1}{1});
    dAnc = max(dAnc, abs(sinrD(pair{1}{2}, iNom) - E2a.KK{j}.viab.SINR_edge_p5));
end
chk('6.3 fila T=1584 vs ffr_results.mat', dAnc, 0);
assert(dAnc == 0, '6.3: el ancla T=1584 NO reproduce ffr_results.mat (%.3e dB).', dAnc);
fprintf('  VERIFICACION OK (ancla exacta, colapso por Delta confirmado)\n');

% --- dibujo ---
COLV   = [0.30 0.30 0.30];   % eje IZQUIERDO (satelites visibles), gris oscuro
COLR   = lines(6);           % MISMOS colores que 6.4/6.6/6.9 para reuse1/D=3/D=4
COLREF = [0.55 0.55 0.55];   % gris de las dos referencias de cfg.viab
COLNOM = [0.45 0.45 0.45];   % gris del marcador del punto nominal

f = figure('Color','w','Visible','off','Position',[100 100 640 440]);
ax = axes(f); hold(ax,'on');

% Eje DERECHO primero: es el que lleva el argumento (KPI) y las dos referencias.
% Cada esquema se distingue por COLOR y ademas por estilo de linea y marcador, para
% que la figura siga siendo legible impresa en escala de grises.
yyaxis(ax,'right');
hR1 = plot(ax, DS.Tlist, sinrD(iR1,:),   '-o', 'Color',COLR(1,:), 'LineWidth',1.6, ...
           'MarkerFaceColor',COLR(1,:), 'MarkerSize',5);
hD3 = plot(ax, DS.Tlist, sinrD(i3(1),:), '--^','Color',COLR(2,:), 'LineWidth',1.6, ...
           'MarkerFaceColor',COLR(2,:), 'MarkerSize',5);
hD4 = plot(ax, DS.Tlist, sinrD(i4(1),:), '-.d','Color',COLR(3,:), 'LineWidth',1.6, ...
           'MarkerFaceColor',COLR(3,:), 'MarkerSize',5);
ylabel(ax, 'SINR_{edge,p5} [dB]');
ax.YAxis(2).Color = 'k';     % tres curvas de colores: el eje va en negro, neutro
ylim(ax, [floor(min([sinrD(iR1,:) gfl6]) - 2.5), ceil(max(sinrD(i4(1),:)) + 6)]);
% Las dos etiquetas van al extremo IZQUIERDO: a la derecha, la de gamma_min chocaba
% con la linea vertical del punto nominal (T=1584). En el extremo izquierdo las dos
% franjas (justo encima de gamma_0 y justo debajo de gamma_min) estan vacias, porque
% a T=66 las tres curvas caen entre -11.6 y -4.6 dB.
yline(ax, gth6, '--', sprintf('\\gamma_0 = %g dB (umbral)', gth6), ...
      'Color',COLREF, 'LineWidth',1.3, 'FontSize',11, 'LabelHorizontalAlignment','left', ...
      'LabelVerticalAlignment','top', 'HandleVisibility','off');
yline(ax, gfl6, '--', sprintf('\\gamma_{min} = %.1f dB (suelo)', gfl6), ...
      'Color',COLREF, 'LineWidth',1.3, 'FontSize',11, 'LabelHorizontalAlignment','left', ...
      'LabelVerticalAlignment','bottom', 'HandleVisibility','off');

% Eje IZQUIERDO. LA ETIQUETA DICE "satelites visibles", NO "interferentes co-canal"
% como decia la figura de FASE A. Lo que este barrido guarda es diag.nVis_mean, los
% satelites de la constelacion sobre la mascara de elevacion, que NO son los
% interferentes co-canal efectivos. La diferencia se ve en el primer punto: a T=66 la
% media vale 0.84 -- MENOR QUE UNO -- porque la cobertura geometrica es 0.774 y el
% promedio incluye los instantes SIN ningun satelite a la vista. Con la etiqueta
% antigua ese 0.84 es ilegible.
yyaxis(ax,'left');
hV = plot(ax, DS.Tlist, nVis, '-s', 'Color',COLV, 'LineWidth',1.6, ...
          'MarkerFaceColor',COLV, 'MarkerSize',5);
ylabel(ax, 'Satelites visibles (media) [-]');
ax.YAxis(1).Color = COLV;
ylim(ax, [0, max(nVis)*1.15]);

% Eje X logaritmico: el eje recorre casi dos decadas (66 -> 4000) y en lineal los
% cuatro primeros puntos se amontonan contra el origen. Marcas SOLO en los siete
% valores barridos, etiquetadas con el valor de T.
set(ax, 'XScale','log', 'XTick', DS.Tlist, 'XMinorTick','off', ...
        'XTickLabel', arrayfun(@(t) sprintf('%d',t), DS.Tlist, 'UniformOutput',false));
xlim(ax, [DS.Tlist(1)*0.85, DS.Tlist(end)*1.15]);
xlabel(ax, 'satelites de la constelacion T [-]');

% Rejilla SOLO vertical: con dos escalas verticales distintas una rejilla horizontal
% pertenece a UNA de las dos y se lee mal. Las dos referencias de cfg.viab hacen de
% unicas reglas horizontales, y son las que importan.
% XMinorGrid OFF explicitamente: en un eje logaritmico MATLAB dibuja la rejilla menor
% por su cuenta al activar XGrid, y aparecian lineas verticales SIN etiqueta en 200,
% 400, 800... justo lo contrario de "marcas solo en los siete valores barridos".
ax.XGrid = 'on';  ax.YGrid = 'off';  ax.XMinorGrid = 'off';
yyaxis(ax,'left');  ax.YGrid = 'off';
box(ax,'on');

% PUNTO DE OPERACION NOMINAL: es lo que ata esta figura al resto del capitulo. La
% Sec. 6.1 trabaja en Walker 53:1584/72/1, y esa fila del barrido es ademas el ancla
% verificada arriba contra ffr_results.mat (max|dif| = 0).
xline(ax, 1584, ':', 'T = 1584 (nominal, Sec. 6.1)', ...
      'Color',COLNOM, 'LineWidth',1.4, 'FontSize',11, 'LabelOrientation','horizontal', ...
      'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','top', ...
      'HandleVisibility','off');

% Leyenda con handles EXPLICITOS: con yyaxis, legend() sin handles recoge solo el
% lado activo y se dejaria fuera la curva del eje izquierdo.
legend(ax, [hV hR1 hD3 hD4], ...
    {'satelites visibles (eje izq.)','reuse1','\Delta = 3','\Delta = 4'}, ...
    'Location','northwest');
bigfont(f); save_fig(f, STAGE, 'fig_6_03_visibles_y_kpi_vs_densidad_orbital'); close(f);
promote('fig_6_03_visibles_y_kpi_vs_densidad_orbital');

%% ---- ARCHIVO: la 6.3 ANTIGUA (FASE A) -- NO es figura del capitulo 6 ----
%  Se conserva y se SIGUE REGENERANDO porque es la figura correcta del PRECEDENTE DE
%  DIMENSIONADO (informe Sec. 7.3): mismo eje T, pero rejilla de 6 km (M=137), Nt=61,
%  dos esquemas y sin veredictos. Se archiva con prefijo ARCHIVO_ -- el mismo criterio
%  con que se conservan las tandas superadas -- para que no se
%  use por error como figura 6.3. Su etiqueta de eje izquierdo ("interferentes
%  co-canal") es la del estudio de FASE A y se deja tal cual: alli T4.nCo si cuenta
%  interferentes co-canal, no satelites visibles, luego la etiqueta es correcta EN SU
%  FIGURA. Es justamente por eso que no se podia reutilizar el nombre.
FA = load(fullfile(PROJROOT,'convergence_studyA.mat'), 'ST','T4');
pub_T   = [66 300 600 1000 1584 2500 4000];
pub_nCo = [0.10 1.45 3.88 6.80 11.56 18.74 30.60];
pub_p5  = [-11.53 -9.98 -6.47 -5.45 -5.04 -4.74 -4.57];

fprintf('\n--- ARCHIVO (FASE A, informe Sec. 7.3) ---\n');
fprintf('  T_list = %s\n', mat2str(FA.ST.T_list));
fprintf('  nCo    = %s\n', mat2str(round(FA.T4.nCo,2)));
fprintf('  p5 r1  = %s\n', mat2str(round(FA.T4.p5(1,:),2)));
assert(isequal(FA.ST.T_list, pub_T), 'ARCHIVO: el eje T no coincide con el publicado.');
dA1 = max(abs(FA.T4.nCo     - pub_nCo)); chk('ARCH n_cocanal vs publicado', dA1, 0.006);
dA2 = max(abs(FA.T4.p5(1,:) - pub_p5 )); chk('ARCH SINR_p5 vs publicado',   dA2, 0.006);
assert(dA1 <= 0.006 && dA2 <= 0.006, 'ARCHIVO no reproduce la tabla de FASE A.');

f = figure('Color','w','Visible','off');
yyaxis left
plot(FA.ST.T_list(FA.T4.ok), FA.T4.nCo(FA.T4.ok), 'o-','LineWidth',1.6);
ylabel('interferentes co-canal visibles [-]');
yyaxis right
plot(FA.ST.T_list(FA.T4.ok), FA.T4.p5(1,FA.T4.ok), 's--','LineWidth',1.4);
ylabel('SINR_{edge,p5} reuse1 [dB]');
grid on; box on; xlabel('satelites de la constelacion T [-]');
bigfont(f); save_fig(f, STAGE, 'ARCHIVO_faseA_interferentes_y_kpi_vs_T_sec7_3'); close(f);
promote('ARCHIVO_faseA_interferentes_y_kpi_vs_T_sec7_3');

%% ==== TABLA DEDUCIDA (Sec. 6.3): geometria co-canal vs ancho de haz ======
%  UNICO bloque del capitulo que NO es CATEGORIA A: no lee ningun .mat de resultados
%  y tampoco simula. Es CALCULO DIRECTO con las funciones del motor --
%    s                  <- build_beam_layout (BL.spacing_deg)
%    distancia de reuso <- ffr_coloring (COL.dmin_s), MEDIDA sobre el layout
%    ganancia           <- beam_gain_dB, la misma que usa el motor
%  -- luego es reproducible sin depender de ninguna tanda. EMITE TABLA POR CONSOLA;
%  no genera figura.
%
%  QUE MUESTRA: que la geometria de la interferencia es INVARIANTE al ancho de haz.
%  Al estrechar el haz, celdas y separacion encogen JUNTAS, asi que la distancia al
%  primer interferente co-canal medida EN ANCHOS DE HAZ no se mueve, y con ella
%  tampoco la ganancia del patron a esa distancia.
%
%  OJO AL ENUNCIADO (es lo que esta tabla precisa): el primer interferente co-canal
%  NO esta a la misma distancia en los tres esquemas. En reuso pleno es el haz
%  ADYACENTE (a s); con reuso Delta esta a la DISTANCIA DE REUSO (sqrt(3)*s con
%  Delta=3, 2*s con Delta=4) y cae bastante mas abajo. Los -10.669 dB del hallazgo de
%  auto-similitud son por tanto los del ADYACENTE, es decir los del reuso pleno.
fprintf('\n--- TABLA DEDUCIDA (Sec. 6.3): geometria co-canal vs ancho de haz ---\n');

TH6  = [3.0 2.2 1.5 1.0 0.7 0.5];        % eje de E3 (run_e3_beamdensity)
cfgT = config_default();
cfgT.beams.nRings = 5;  cfgT.beams.nBeams = [];   % 91 haces: cluster de sobra para
                                                  % que el par co-canal mas proximo
                                                  % de cada Delta exista en el layout

% (i) DISTANCIA DE REUSO: la MIDE ffr_coloring sobre el layout (COL.dmin_s); no se
%     escribe a mano. Si dejara de valer sqrt(3) y 2 seria un RESULTADO y no un
%     detalle, asi que se comprueba con assert en vez de asumirse.
nK   = numel(TH6);
sdeg = nan(1,nK);  k3 = nan(1,nK);  k4 = nan(1,nK);
for k = 1:nK
    c = cfgT;  c.radio.beamwidth3dB_deg = TH6(k);
    BLk = build_beam_layout(c);
    [~, C3] = ffr_coloring(BLk, 3, false);
    [~, C4] = ffr_coloring(BLk, 4, false);
    sdeg(k) = BLk.spacing_deg;  k3(k) = C3.dmin_s;  k4(k) = C4.dmin_s;
end
fprintf('  distancia de reuso MEDIDA por ffr_coloring (en unidades de s):\n');
fprintf('    Delta=3 -> %.10f  (sqrt(3) = %.10f)   max|dif| = %.2e\n', ...
    mean(k3), sqrt(3), max(abs(k3 - sqrt(3))));
fprintf('    Delta=4 -> %.10f  (2        = %.10f)   max|dif| = %.2e\n', ...
    mean(k4), 2, max(abs(k4 - 2)));
assert(max(abs(k3-sqrt(3))) < 1e-9 && max(abs(k4-2)) < 1e-9, ...
    ['Sec 6.3: la distancia de reuso MEDIDA no es sqrt(3)*s / 2*s (D=3: %.6f, ' ...
     'D=4: %.6f). Es un RESULTADO, no un detalle: revisar antes de publicar.'], ...
    mean(k3), mean(k4));

% (ii) Ganancia relativa del patron a cada distancia. SIN suelo de lobulos: la tabla
%      describe el PATRON. Se comprueba abajo que el suelo no morderia de todos modos.
G1 = arrayfun(@(k) beam_gain_dB(        sdeg(k), TH6(k)), 1:nK);   % adyacente,     s
G3 = arrayfun(@(k) beam_gain_dB(k3(k).*sdeg(k), TH6(k)), 1:nK);    % Delta=3, sqrt3*s
G4 = arrayfun(@(k) beam_gain_dB(k4(k).*sdeg(k), TH6(k)), 1:nK);    % Delta=4,     2*s

% (iii) LIMITE de angulo pequeno (theta -> 0), ANALITICO. Esta FUERA del eje barrido.
%       Con s_rad -> sqrt(3)*sin(theta/2), el argumento u = 1.61634*sin(d)/sin(theta/2)
%       tiende a un valor FIJO por distancia y las seis columnas se hacen constantes:
%         d = s          -> u = 1.61634*sqrt(3)   ;  s/theta       -> sqrt(3)/2
%         d = sqrt(3)*s  -> u = 1.61634*3         ;  sqrt(3)*s/th  -> 3/2
%         d = 2*s        -> u = 1.61634*2*sqrt(3) ;  2*s/theta     -> sqrt(3)
uL = 1.61634 * [sqrt(3), 3, 2*sqrt(3)];
GL = 20*log10(abs(2*besselj(1,uL)./uL));
rL = [sqrt(3)/2, 3/2, sqrt(3)];

fprintf('\n%10s %10s | %12s %10s | %12s %10s | %12s %10s\n', ...
    'theta_3dB','s','s/theta','G(s)','sqrt3*s/th','G(sqrt3*s)','2*s/theta','G(2*s)');
fprintf('%10s %10s | %12s %10s | %12s %10s | %12s %10s\n', ...
    '[deg]','[deg]','[anchos]','[dB]','[anchos]','[dB]','[anchos]','[dB]');
for k = 1:nK
    fprintf('%10.2f %10.4f | %12.4f %10.3f | %12.4f %10.3f | %12.4f %10.3f\n', ...
        TH6(k), sdeg(k), sdeg(k)/TH6(k), G1(k), ...
        k3(k)*sdeg(k)/TH6(k), G3(k), k4(k)*sdeg(k)/TH6(k), G4(k));
end
fprintf('%10s %10s | %12.4f %10.3f | %12.4f %10.3f | %12.4f %10.3f\n', ...
    'lim th->0','   --   ', rL(1), GL(1), rL(2), GL(2), rL(3), GL(3));
fprintf('  (ultima fila = LIMITE ANALITICO de angulo pequeno, FUERA del eje barrido)\n');
% PIE DE TABLA (hay que llevarlo a la memoria junto con la tabla): todas las
% distancias son ENTRE CENTROS DE CELDA. El usuario de BORDE no esta en el centro:
% esta desplazado hasta el circunradio s/sqrt(3) = 0.5 anchos de haz, y desde ahi ve
% OTRAS distancias a sus co-canal -- y no las mismas en los dos esquemas, porque la
% sub-reticula de Delta=3 esta girada 30 deg (sus co-canal caen en la direccion del
% VERTICE) y la de Delta=4 no (caen en la del ADYACENTE). Por eso esta tabla NO
% explica por si sola el hueco de SINR_edge_p5 entre Delta=3 y Delta=4; ese analisis
% es del capitulo 8.
fprintf('  PIE: las distancias son ENTRE CENTROS DE CELDA. El usuario de BORDE esta\n');
fprintf('       desplazado hasta el circunradio (s/sqrt(3) = 0.5 anchos) y ve otras\n');
fprintf('       distancias a sus co-canal, distintas ademas en cada esquema (cap. 8).\n');

% --- comprobacion 1: s/theta tiende a sqrt(3)/2; donde deja de ser constante ---
%     La cifra se BUSCA, no se deduce del orden de magnitud: dos valores separados
%     9.6e-05 pueden diferir ya en la 4a decimal si caen a distinto lado del redondeo.
r1 = sdeg./TH6;
dfirst = NaN;
for d = 1:12
    if round(r1(1),d) ~= round(r1(end),d), dfirst = d; break; end
end
fprintf('\n  [1] s/theta: %.10f (th=3.0) -> %.10f (th=0.5) | limite sqrt(3)/2 = %.10f\n', ...
    r1(1), r1(end), sqrt(3)/2);
fprintf('      recorrido = %.3e | constante en %d decimales; la PRIMERA cifra que se\n', ...
    abs(r1(1)-r1(end)), dfirst-1);
fprintf('      mueve es la %d.a: %.4f (th=3.0) vs %.4f (th=0.5)\n', dfirst, r1(1), r1(end));

% --- comprobacion 2: la cadena es la correcta (ancho nominal -> memoria) ---
cN = cfgT;  cN.radio.beamwidth3dB_deg = 2.08;
BLN  = build_beam_layout(cN);
proj = cN.constellations(1).h * tand(BLN.spacing_deg);
fprintf('\n  [2] ancho nominal 2.08 deg -> s = %.6f deg (memoria 1.8012) | %g*tand(s) = %.4f km\n', ...
    BLN.spacing_deg, cN.constellations(1).h, proj);
d63s = max(abs(BLN.spacing_deg - 1.8012), abs(proj - 17.30)/100);
chk('6.3 cadena s / paso proyectado', d63s, 5e-5);
assert(abs(BLN.spacing_deg - 1.8012) < 5e-5 && abs(proj - 17.30) < 5e-3, ...
    'Sec 6.3: la cadena no reproduce s = 1.8012 deg / 17.30 km (memoria).');

% --- comprobacion 3: proteccion del reuso, MEDIDA (no deducida) ---
p3 = G1 - G3;  p4 = G1 - G4;
fprintf('\n  [3] proteccion del reuso = G(adyacente) - G(primer co-canal):\n');
fprintf('      Delta=3 : %.3f .. %.3f dB (media %.3f) | limite th->0 %.3f dB\n', ...
    min(p3), max(p3), mean(p3), GL(1)-GL(2));
fprintf('      Delta=4 : %.3f .. %.3f dB (media %.3f) | limite th->0 %.3f dB\n', ...
    min(p4), max(p4), mean(p4), GL(1)-GL(3));

% --- comprobacion 4: el suelo de lobulos NO interviene en esta tabla ---
fprintf('\n  [4] suelo cfg.radio.sidelobe_floor_dB = %.1f dB | minimo de la tabla = %.3f dB\n', ...
    cfgT.radio.sidelobe_floor_dB, min([G1 G3 G4]));
assert(min([G1 G3 G4]) > cfgT.radio.sidelobe_floor_dB, ...
    'Sec 6.3: el suelo de lobulos morderia en esta tabla; ya no es patron puro.');
fprintf('      el suelo NO muerde: la tabla es patron puro\n');

% --- comprobacion 5: control cruzado del limite con la funcion del motor ---
thTiny = 1e-4;  sTiny = rad2deg(sqrt(3)*sin(deg2rad(thTiny)/2));
Gt = [beam_gain_dB(sTiny,thTiny), beam_gain_dB(sqrt(3)*sTiny,thTiny), ...
      beam_gain_dB(2*sTiny,thTiny)];
dLim = max(abs(Gt - GL));
chk('6.3 limite analitico vs beam_gain_dB', dLim, 1e-6);
assert(dLim < 1e-6, 'Sec 6.3: el limite analitico no coincide con beam_gain_dB.');
fprintf('  VERIFICACION OK\n');

%% ============ 6.4 y 6.5  E3: densidad de haces ===========================
E3 = load(fullfile(PROJROOT,'e3_results.mat'), 'resA','resB','dens','biasedA','cases','cfg');
nC3 = numel(E3.cases);  nT3 = numel(E3.dens);
nm3 = cellfun(@(c) c.name, E3.cases, 'UniformOutput', false);
col3 = lines(nC3);
gth = E3.cfg.viab.gamma_th_dB;  gfl = E3.cfg.viab.gamma_floor_dB;
sinrB = @(c) arrayfun(@(i) E3.resB{i}.K{c}.viab.SINR_edge_p5, 1:nT3);
sinrA = @(c) arrayfun(@(i) E3.resA{i}.K{c}.viab.SINR_edge_p5, 1:nT3);

% Tabla publicada de E3 variante B, en dens = 1.854 / 7.424 / 16.706 / 66.830
iSel   = [1 3 4 6];
pub64  = [-4.91 -4.93 -5.21 -5.76;    % reuse1
           3.93   3.97  3.25  2.04;   % ffr(3,a*)  -> indice 4
           7.56   7.41  6.58  4.56];  % ffr(4,a*)  -> indice 5
got64  = [sinrB(1); sinrB(4); sinrB(5)];
fprintf('\n--- 6.4 / 6.5 (E3) ---\n');
fprintf('  dens     = %s\n', mat2str(round(E3.dens,3)));
fprintf('  reuse1 B = %s\n', mat2str(round(sinrB(1),2)));
fprintf('  ffr(3) B = %s\n', mat2str(round(sinrB(4),2)));
fprintf('  ffr(4) B = %s\n', mat2str(round(sinrB(5),2)));
d64 = max(max(abs(got64(:,iSel) - pub64)));
chk('6.4 SINR_edge_p5 (var. B) vs publ.', d64, 0.006);
assert(d64 <= 0.006, '6.4 no reproduce la tabla de E3.');
% veredictos publicados: reuse1 marginal en TODA la banda; ffr(3)/ffr(4) viable
vB1 = arrayfun(@(i) string(E3.resB{i}.K{1}.viab.verdict), 1:nT3);
vB4 = arrayfun(@(i) string(E3.resB{i}.K{4}.viab.verdict), 1:nT3);
fprintf('  veredicto reuse1 (B) = %s\n', strjoin(cellstr(vB1),' '));
fprintf('  veredicto ffr(3) (B) = %s\n', strjoin(cellstr(vB4),' '));
assert(all(vB1 == "marginal"), '6.5: reuse1 deberia ser marginal en todo el eje.');
assert(all(vB4 == "viable"),   '6.5: ffr(3) deberia ser viable en todo el eje.');
fprintf('  veredictos OK\n');

% --- 6.4 ---
%  Leyenda COMPARTIDA fuera de los ejes (tile 'south'): con 6 esquemas, dentro del
%  area de dibujo tapaba curvas (en la 6.6 llegaba a ocultar la de reuse1).
f = figure('Color','w','Visible','off','Position',WL);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
axL = gobjects(1,2);
for v = 1:2
    ax = nexttile(tl); axL(v) = ax; hold on; grid on; box on;
    for c = 1:nC3
        if v == 1, yy = sinrA(c); else, yy = sinrB(c); end
        plot(E3.dens, yy, '-o', 'Color', col3(c,:), 'LineWidth',1.4, ...
             'MarkerFaceColor', col3(c,:), 'MarkerSize',3.5, 'DisplayName', nm3{c});
    end
    yline(gth,'k--','LineWidth',1.2,'HandleVisibility','off');
    yline(gfl,'r--','LineWidth',1.2,'HandleVisibility','off');
    if v == 1
        for i = find(E3.biasedA)
            xline(E3.dens(i), ':', 'Color',[.6 .6 .6], 'HandleVisibility','off');
        end
    end
    set(gca,'XScale','log');
    xlabel('densidad de haces [haces/1000 km^2]'); ylabel('SINR_{edge,p5} [dB]');
    panel(ax, sprintf('(%c)', 'a'+v-1));
end
lgd = legend(axL(1), nm3, 'Orientation','horizontal','NumColumns',3);
lgd.Layout.Tile = 'south';
bigfont(f); save_fig(f, STAGE, 'fig_6_04_sinr_borde_vs_densidad'); close(f);
promote('fig_6_04_sinr_borde_vs_densidad');

% --- 6.5 mapa de veredictos ---
lev = containers.Map({'inviable','marginal','viable','sin-datos'},{1,2,3,0});
f = figure('Color','w','Visible','off','Position',W2);
for v = 1:2
    ax = subplot(1,2,v);
    if v==1, res = E3.resA; else, res = E3.resB; end
    Mv = zeros(nC3, nT3);
    for c = 1:nC3, for i = 1:nT3, Mv(c,i) = lev(res{i}.K{c}.viab.verdict); end, end
    imagesc(Mv,[0 3]);
    colormap([0.85 0.85 0.85; 0.80 0.20 0.20; 0.95 0.75 0.20; 0.20 0.65 0.30]);
    set(gca,'YTick',1:nC3,'YTickLabel',nm3);
    set(gca,'XTick',1:nT3,'XTickLabel',arrayfun(@(d) sprintf('%.2f',d), E3.dens,'UniformOutput',false));
    set(gca,'XTickLabelRotation',45);
    xlabel('densidad de haces [haces/1000 km^2]');
    if v==1
        % 'SESG' VERTICAL dentro de la fila 1: en horizontal y a 11 pt las tres
        % columnas sesgadas (34.10 / 66.83 y la contigua) se solapaban entre si.
        for i = find(E3.biasedA)
            text(i, 1, 'SESG', 'Rotation',90, 'HorizontalAlignment','center', ...
                 'VerticalAlignment','middle', 'FontSize',11, 'FontWeight','bold', ...
                 'Color','w');
        end
    end
    panelOut(ax, sprintf('(%c)', 'a'+v-1));
end
bigfont(f); save_fig(f, STAGE, 'fig_6_05_mapa_veredictos_densidad'); close(f);
promote('fig_6_05_mapa_veredictos_densidad');

%% ============ 6.6 y 6.7  E3b: mascara de elevacion =======================
EB = load(fullfile(PROJROOT,'minelev_results.mat'), 'R','minElevList','densList','nm','cases','cfg');
nD = numel(EB.densList);  nE = numel(EB.minElevList);  nCb = numel(EB.nm);
colb = lines(nCb);
gthb = EB.cfg.viab.gamma_th_dB;  gflb = EB.cfg.viab.gamma_floor_dB;
sB = @(d,c) arrayfun(@(i) EB.R{d}{i}.K{c}.viab.SINR_edge_p5, 1:nE);

fprintf('\n--- 6.6 / 6.7 (E3b) ---\n');
fprintf('  densidades: %s (T=%d) | %s (T=%d)\n', EB.densList(1).name, EB.densList(1).T, ...
        EB.densList(2).name, EB.densList(2).T);
fprintf('  minElev   = %s\n', mat2str(EB.minElevList));
% publicado (T=66, densList(2)): minElev 45 40 35 30 25 ... 10
iP = [1 2 3 4 5 8];
pub66 = [-5.66 -7.41 -8.37  -9.91 -11.74 -14.93;   % reuse1
          2.29 -0.73 -2.17  -4.20  -6.38  -9.90;   % reuseD(3)
          5.45  1.82  0.08  -2.39  -4.81  -8.56];  % reuseD(4)
got66 = [sB(2,1); sB(2,2); sB(2,3)];
fprintf('  T=66 reuse1    = %s\n', mat2str(round(sB(2,1),2)));
fprintf('  T=66 reuseD(3) = %s\n', mat2str(round(sB(2,2),2)));
d66 = max(max(abs(got66(:,iP) - pub66)));
chk('6.6 SINR_edge_p5 (T=66) vs publ.', d66, 0.006);
assert(d66 <= 0.006, '6.6 no reproduce la tabla de E3b.');
% el barrido a T=1584 debe ser DEGENERADO (recorrido 0.000 dB en las 8 mascaras)
% DEGENERACION del barrido denso. Lo publicado dice "recorrido = 0.000 dB en las 8
% mascaras (filas BIT-IDENTICAS)". Medido, el recorrido NO es exactamente cero:
% es del orden de 1e-6 dB. La cifra publicada es correcta a 3 decimales, pero
% "bit-identicas" es una exageracion. Se mide, se imprime y se reporta.
recAll = arrayfun(@(c) max(sB(1,c)) - min(sB(1,c)), 1:nCb);
fprintf('  recorrido a T=1584 por esquema [dB]:\n');
for c = 1:nCb
    fprintf('      %-14s %.3e\n', EB.nm{c}, recAll(c));
end
fprintf('  recorrido MAXIMO = %.3e dB (publicado "0.000 dB, filas bit-identicas")\n', max(recAll));
assert(max(recAll) < 1e-5, ...
    '6.6: el barrido denso no es degenerado (recorrido %.3e dB).', max(recAll));
fprintf('  degeneracion OK a 3 decimales (pero NO bit-identica: ver informe)\n');

f = figure('Color','w','Visible','off','Position',WL);
tl = tiledlayout(f,1,nD,'TileSpacing','compact','Padding','compact');
axL = gobjects(1,nD);
for d = 1:nD
    ax = nexttile(tl); axL(d) = ax; hold on; grid on; box on;
    for c = 1:nCb
        plot(EB.minElevList, sB(d,c), '-o', 'Color',colb(c,:), 'LineWidth',1.4, ...
             'MarkerFaceColor',colb(c,:), 'MarkerSize',3.5, 'DisplayName', EB.nm{c});
    end
    yline(gthb,'k--','LineWidth',1.2,'HandleVisibility','off');
    yline(gflb,'r--','LineWidth',1.2,'HandleVisibility','off');
    set(gca,'XDir','reverse');
    xlabel('mascara de elevacion [deg]'); ylabel('SINR_{edge,p5} [dB]');
    panel(ax, sprintf('(%c)', 'a'+d-1));
end
lgd = legend(axL(1), EB.nm, 'Orientation','horizontal','NumColumns',3);
lgd.Layout.Tile = 'south';
bigfont(f); save_fig(f, STAGE, 'fig_6_06_sinr_borde_vs_mascara'); close(f);
promote('fig_6_06_sinr_borde_vs_mascara');

f = figure('Color','w','Visible','off','Position',W2);
for d = 1:nD
    ax = subplot(1,nD,d);
    Mv = zeros(nCb, nE);
    for c = 1:nCb, for i = 1:nE, Mv(c,i) = lev(EB.R{d}{i}.K{c}.viab.verdict); end, end
    imagesc(Mv,[0 3]);
    colormap([0.85 0.85 0.85; 0.80 0.20 0.20; 0.95 0.75 0.20; 0.20 0.65 0.30]);
    set(gca,'YTick',1:nCb,'YTickLabel',EB.nm);
    set(gca,'XTick',1:nE,'XTickLabel',arrayfun(@(x) sprintf('%d',x), EB.minElevList,'UniformOutput',false));
    xlabel('mascara de elevacion [deg]');
    % Abreviatura del veredicto en cada celda (como el original): sin ella el mapa
    % solo se lee por color y pierde informacion al imprimirse en escala de grises.
    for c = 1:nCb
        for i = 1:nE
            text(i, c, EB.R{d}{i}.K{c}.viab.verdict(1:3), 'HorizontalAlignment','center', ...
                 'FontSize',11, 'Color','w', 'FontWeight','bold');
        end
    end
    panelOut(ax, sprintf('(%c)', 'a'+d-1));
end
bigfont(f); save_fig(f, STAGE, 'fig_6_07_mapa_veredictos_mascara'); close(f);
promote('fig_6_07_mapa_veredictos_mascara');

%% ====== ANEXO: 6.8 antigua (KPI vs numero de operadores) ================
%  NO es figura del cuerpo: la sustituyo la 6.9. Se REGENERA igualmente para que el
%  catalogo siga cumpliendo lo que promete -- que TODAS sus figuras salen de los .mat
%  y ninguna es un PNG huerfano. Es el mismo criterio con el que se conserva
%  ARCHIVO_faseA_... mas arriba. Bloque recuperado tal cual, con sus verificaciones.
E8 = load(fullfile(PROJROOT,'interconstellation_results.mat'), 'res','nOpsList','nm','cases','cfg');
nQ8 = numel(E8.nOpsList);  nC8 = numel(E8.nm);
col8 = lines(nC8);
p1_8 = @(c) arrayfun(@(q) E8.res{q}.K{c}.viab.SINR_edge_p1, 1:nQ8);
p5_8 = @(c) arrayfun(@(q) E8.res{q}.K{c}.viab.SINR_edge_p5, 1:nQ8);

fprintf('\n--- ANEXO 6.8 (E4, figura de anexo) ---\n');
pub8P1 = [-5.407 2.834 6.042 2.834 6.042 2.834];
pub8P5 = [-4.950 3.865 7.226 3.865 7.226 3.865];
got8P1 = arrayfun(@(c) E8.res{1}.K{c}.viab.SINR_edge_p1, 1:nC8);
got8P5 = arrayfun(@(c) E8.res{1}.K{c}.viab.SINR_edge_p5, 1:nC8);
d68a = max(abs(got8P1 - pub8P1)); chk('6.8 SINR_edge_p1 vs publicado', d68a, 0.0006);
d68b = max(abs(got8P5 - pub8P5)); chk('6.8 SINR_edge_p5 vs publicado', d68b, 0.0006);
dInv = max([ max(arrayfun(@(c) max(abs(diff(p1_8(c)))), 1:nC8)), ...
             max(arrayfun(@(c) max(abs(diff(p5_8(c)))), 1:nC8)) ]);
chk('6.8 invariancia 1->2 operadores', dInv, 5e-4);
assert(d68a <= 0.0006 && d68b <= 0.0006, '6.8 no reproduce la tabla de E4.');
assert(dInv < 5e-4, '6.8: la invariancia 1->2 operadores no se cumple (%.6f dB).', dInv);

f = figure('Color','w','Visible','off','Position',WL);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
axL = gobjects(1,2);
for v = 1:2
    ax = nexttile(tl); axL(v) = ax; hold on; grid on; box on;
    for c = 1:nC8
        if v==1, yy = p1_8(c); else, yy = p5_8(c); end
        plot(E8.nOpsList, yy, '-o','Color',col8(c,:),'LineWidth',1.4, ...
             'MarkerFaceColor',col8(c,:),'MarkerSize',4,'DisplayName',E8.nm{c});
    end
    yline(E8.cfg.viab.gamma_th_dB,'k--','HandleVisibility','off');
    yline(E8.cfg.viab.gamma_floor_dB,'r--','HandleVisibility','off');
    xlabel('numero de operadores co-canal [-]');
    if v==1, ylabel('SINR_{edge,p1} [dB]'); else, ylabel('SINR_{edge,p5} [dB]'); end
    set(gca,'XTick',E8.nOpsList);
    panel(ax, sprintf('(%c)', 'a'+v-1));
end
lgd = legend(axL(1), E8.nm, 'Orientation','horizontal','NumColumns',3);
lgd.Layout.Tile = 'south';
bigfont(f); save_fig(f, STAGE, 'fig_6_08_kpi_vs_operadores'); close(f);
promote('fig_6_08_kpi_vs_operadores');
fprintf('  regenerada como material de ANEXO (la figura del cuerpo es la 6.9)\n');

%% ====== 6.9  E4: eventos in-line co-canal y veredicto vs operadores ======
%  SUSTITUYE a la antigua fig_6_08_kpi_vs_operadores (SINR_edge_p1 y p5 frente al
%  numero de operadores). Aquella gastaba dos paneles para dos puntos y tres rectas
%  horizontales: su contenido completo cabe en una tabla de dos filas y no anadia
%  nada a la prosa. El panel (a) nuevo aporta lo que la tabla no puede dar -- la
%  DISTRIBUCION de la separacion angular al satelite ajeno mejor alineado --, que es
%  precisamente la razon de que el veredicto no se mueva al anadir un operador.
%  fig_6_08_kpi_vs_operadores.png NO se borra: pasa a material de ANEXO.
%
%  CATEGORIA A: todo se LEE de interconstellation_results.mat. No se llama a
%  compute_sinr_ffr / compute_interference / compute_kpis / run_sweep_points.
E4 = load(fullfile(PROJROOT,'interconstellation_results.mat'), ...
          'INL','res','nOpsList','nm','cases','cfg');
nQ = numel(E4.nOpsList);  nC4 = numel(E4.nm);
col4  = lines(nC4);                       % MISMOS colores que 6.4 y 6.6
FIG69 = 'fig_6_09_operadores_cocanal';    % el nombre, en UN SOLO sitio
psiTh = [1 2 5 10];                       % deg, los umbrales de eventos de E4.
                                          % No se leen del .mat a proposito: la
                                          % comprobacion (3) contra INL.P_event
                                          % falla si no son los mismos.
p1 = @(c) arrayfun(@(q) E4.res{q}.K{c}.viab.SINR_edge_p1, 1:nQ);
p5 = @(c) arrayfun(@(q) E4.res{q}.K{c}.viab.SINR_edge_p5, 1:nQ);

fprintf('\n--- 6.9 (E4: eventos in-line + veredicto) ---\n');
fprintf('  nOps = %s\n', mat2str(E4.nOpsList));

% (1) SINR_edge_p5 de los 6 esquemas con 1 operador, contra lo publicado
pubP5 = [-4.950 3.865 7.226 3.865 7.226 3.865];
pubP1 = [-5.407 2.834 6.042 2.834 6.042 2.834];
gotP5 = arrayfun(@(c) E4.res{1}.K{c}.viab.SINR_edge_p5, 1:nC4);
gotP1 = arrayfun(@(c) E4.res{1}.K{c}.viab.SINR_edge_p1, 1:nC4);
fprintf('  p5 (1 op) = %s\n', mat2str(round(gotP5,3)));
fprintf('  p1 (1 op) = %s   <- va al TEXTO, ya no a la figura\n', mat2str(round(gotP1,3)));
d69a = max(abs(gotP5 - pubP5));  chk('6.9 (1) SINR_edge_p5 vs publicado', d69a, 0.0006);
assert(d69a <= 0.0006, '6.9: SINR_edge_p5 no reproduce la tabla de E4 (%.6f dB).', d69a);
% Guarda HEREDADA del bloque anterior: el p1 sale de la figura pero sigue yendo al
% texto, asi que se conserva su comprobacion contra lo publicado.
d69a1 = max(abs(gotP1 - pubP1)); chk('6.9 (1b) SINR_edge_p1 vs publicado', d69a1, 0.0006);
assert(d69a1 <= 0.0006, '6.9: SINR_edge_p1 no reproduce la tabla de E4 (%.6f dB).', d69a1);

% (2) invariancia 1 -> 2 operadores, en p5 Y en p1 de los 6 esquemas
d69b = max([ max(arrayfun(@(c) max(abs(diff(p5(c)))), 1:nC4)), ...
             max(arrayfun(@(c) max(abs(diff(p1(c)))), 1:nC4)) ]);
chk('6.9 (2) invariancia 1->2 operadores', d69b, 5e-4);
assert(d69b < 5e-4, '6.9: la invariancia 1->2 operadores no se cumple (%.3e dB).', d69b);

% (3) La CURVA que se dibuja tiene que ser la MISMA poblacion que la tabla de
%     eventos: se recalcula P(psi < X) desde el vector ordenado que se grafica y se
%     contrasta contra INL{2}.P_event. Sin esto, el panel podria dibujar una
%     poblacion distinta de la que sostiene la tabla y nadie se enteraria.
%     Con 1 operador NO hay satelites ajenos (psi_min todo NaN), por eso solo INL{2}.
E   = E4.INL{2};
valid = E.cov & ~isnan(E.psi_min);
x   = sort(E.psi_min(valid));   nV = numel(x);
ycdf = (1:nV)/nV;
pubPev = [0.00027774 0.00222191 0.03110679 0.12234412];
gotPev = arrayfun(@(t) sum(x < t)/nV, psiTh);
fprintf('  P(psi<%g) P(psi<%g) P(psi<%g) P(psi<%g) = %s\n', psiTh(1), psiTh(2), ...
    psiTh(3), psiTh(4), mat2str(round(gotPev,8)));
d69c = max(abs(gotPev - pubPev));  chk('6.9 (3) CDF reproduce INL.P_event', d69c, 1e-6);
assert(d69c <= 1e-6, ...
    '6.9: la curva dibujada NO es la poblacion de la tabla de eventos (%.3e).', d69c);

% (4) la poblacion en si: numero de instantes validos y estadisticos de psi_min
fprintf('  n valid = %d | psi_min medio = %.3f deg | psi_min min = %.3f deg\n', ...
    nV, mean(x), min(x));
d69d = max([abs(nV - 7201), abs(mean(x) - 22.368), abs(min(x) - 0.362)]);
chk('6.9 (4) poblacion de psi_min', d69d, 0.001);
assert(d69d <= 0.001, '6.9: la poblacion de psi_min no coincide (%.4f).', d69d);
fprintf('  VERIFICACION OK (4/4) -> se promueve la figura\n');

% ------------------------------ dibujo ----------------------------------
f = figure('Color','w','Visible','off','Position',WL);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

% --- (a) CDF de psi_min con DOS operadores ---
%  EJE Y LOGARITMICO A PROPOSITO, y es el motivo de que el panel exista: la cola
%  informativa va de 2.8e-4 a 0.12 y en eje lineal queda aplastada contra el suelo.
%  En log se leen las tres decadas donde viven los eventos.
%  SIN LEYENDA, y es deliberado: con 1 operador no hay satelites ajenos (psi_min todo
%  NaN) y no se dibuja ninguna curva para ese caso. El make_e4_figures original
%  construia igualmente la leyenda con las DOS etiquetas, y MATLAB asignaba "1 op" a
%  la unica curva dibujada, que es la de DOS operadores: la leyenda de aquel PNG es
%  FALSA. Aqui hay una sola curva y el caso va en el pie de figura.
axA = nexttile(tl); hold(axA,'on'); grid(axA,'on'); box(axA,'on');
plot(axA, x, ycdf, 'LineWidth',1.8, 'Color',[0.20 0.45 0.70]);
set(axA, 'XScale','log', 'YScale','log');
ylim(axA, [1e-4 1]);
for j = 1:numel(psiTh)
    xline(axA, psiTh(j), ':', sprintf('%g',psiTh(j)), 'Color',[0.45 0.45 0.45], ...
          'LineWidth',1.2, 'FontSize',11, 'LabelVerticalAlignment','bottom', ...
          'LabelHorizontalAlignment','left', 'HandleVisibility','off');
end
xlabel(axA, 'separacion angular al ajeno mejor alineado \psi_{min} [deg]');
ylabel(axA, 'F(\psi_{min}) [-]');
panel(axA, '(a)');

% --- (b) indicador de borde frente al numero de operadores ---
%  Solo SINR_edge_p5. El percentil 1 se retira de la figura y va al texto.
axB = nexttile(tl); hold(axB,'on'); grid(axB,'on'); box(axB,'on');
for c = 1:nC4
    plot(axB, E4.nOpsList, p5(c), '-o', 'Color',col4(c,:), 'LineWidth',1.4, ...
         'MarkerFaceColor',col4(c,:), 'MarkerSize',5, 'DisplayName',E4.nm{c});
end
yline(axB, E4.cfg.viab.gamma_th_dB,    'k--','LineWidth',1.2,'HandleVisibility','off');
yline(axB, E4.cfg.viab.gamma_floor_dB, 'r--','LineWidth',1.2,'HandleVisibility','off');
xlabel(axB, 'numero de operadores co-canal [-]');
ylabel(axB, 'SINR_{edge,p5} [dB]');
set(axB, 'XTick', E4.nOpsList);
panel(axB, '(b)');

% La leyenda cuelga del eje del panel (b): el (a) ya no tiene series con nombre.
lgd = legend(axB, E4.nm, 'Orientation','horizontal','NumColumns',3);
lgd.Layout.Tile = 'south';
bigfont(f); save_fig(f, STAGE, FIG69); close(f);
promote(FIG69);

%% ====== 6.10  Mapa de veredictos UNIFICADO de los cuatro ejes ===========
%  La figura que faltaba: los CUATRO ejes de saturacion del capitulo en una sola
%  imagen, 2x2, cada panel una rejilla de 3 filas (reuse1 / D=3 / D=4) por tantas
%  columnas como puntos tenga el eje. Es el argumento del capitulo hecho imagen: la
%  fila de reuse1 no se pone VERDE en ninguna celda de ningun eje.
%  CATEGORIA A: se LEE de los cuatro .mat; ninguna llamada al motor.
%  Estilo CALCADO de los mapas 6.5 y 6.7: mismo `lev`, misma colormap, imagesc con
%  limites [0 3], abreviatura de 3 letras en blanco y negrita dentro de cada celda y
%  panelOut para la etiqueta, que aqui va FUERA del eje porque las celdas llegan al
%  borde.
fprintf('\n--- 6.10 (mapa de veredictos unificado) ---\n');

DS10 = load(fullfile(PROJROOT,'e3_density_sweep.mat'),          'results','Tlist','schemeCases');
E310 = load(fullfile(PROJROOT,'e3_results.mat'),                'resB','dens','cases');
EB10 = load(fullfile(PROJROOT,'minelev_results.mat'),           'R','minElevList','densList','cases');
E410 = load(fullfile(PROJROOT,'interconstellation_results.mat'),'res','nOpsList','cases');

% R{2} de E3b es la constelacion DISPERSA (T=66), que es la que instrumenta el eje;
% R{1} (densa) da las ocho columnas identicas y no aporta nada a este mapa.
assert(EB10.densList(2).T == 66, ...
    '6.10: R{2} de minelev_results no es la constelacion dispersa T=66 (T=%d).', EB10.densList(2).T);

EJ = struct( ...
  'tag',  {'(a)','(b)','(c)','(d)'}, ...
  'R',    {DS10.results, E310.resB, EB10.R{2}, E410.res}, ...
  'cases',{DS10.schemeCases, E310.cases, EB10.cases, E410.cases}, ...
  'lab',  {arrayfun(@(t) sprintf('%d',t),     DS10.Tlist,        'UniformOutput',false), ...
           arrayfun(@(d) sprintf('%.2f',d),   E310.dens,         'UniformOutput',false), ...
           arrayfun(@(e) sprintf('%d',e),     EB10.minElevList,  'UniformOutput',false), ...
           arrayfun(@(q) sprintf('%d',q),     E410.nOpsList,     'UniformOutput',false)}, ...
  'xlab', {'satelites de la constelacion T [-]', 'densidad de haces [haces/1000 km^2]', ...
           'mascara de elevacion [deg]',         'numero de operadores co-canal [-]'}, ...
  'rot',  {45, 45, 45, 0});

% --- (0) COLAPSO DE ESQUEMAS: es lo que LEGITIMA la figura, asi que va primero ----
%  Los .mat traen 6 esquemas y el mapa dibuja 3 filas. Eso solo vale si dentro de
%  cada Delta el veredicto coincide entre reuso-Delta, FFR estatica y FFR adaptativa.
%  Los grupos se derivan del campo Delta de cada caso -- NO del orden, que aunque hoy
%  es el mismo en los cuatro ficheros no tiene por que serlo.
nMis = 0;  nChk = 0;  repAll = nan(4,3);
for e = 1:4
    Dl = cellfun(@(c) c.Delta, EJ(e).cases);
    g1 = find(Dl == 1);  g3 = find(Dl == 3);  g4 = find(Dl == 4);
    assert(numel(g1) == 1, '6.10: el eje %s no tiene exactamente un esquema reuse1.', EJ(e).tag);
    assert(~isempty(g3) && ~isempty(g4), '6.10: falta algun Delta en el eje %s.', EJ(e).tag);
    for i = 1:numel(EJ(e).R)
        for g = {g3, g4}
            vv = arrayfun(@(c) string(EJ(e).R{i}.K{c}.viab.verdict), g{1});
            nMis = nMis + sum(vv ~= vv(1));
            nChk = nChk + numel(vv) - 1;
        end
    end
    repAll(e,:) = [g1, g3(1), g4(1)];       % representantes de cada fila
    fprintf('  %s grupos por Delta: reuse1=%d | D=3=%s | D=4=%s\n', EJ(e).tag, g1, ...
        mat2str(g3(:).'), mat2str(g4(:).'));
end
chk('6.10 (0) colapso de esquemas por Delta', nMis, 0);
assert(nMis == 0, ...
    ['6.10: hay %d celdas donde los esquemas del MISMO Delta NO comparten veredicto. ' ...
     'La figura de 3 filas no seria legitima: es un RESULTADO, no un detalle.'], nMis);
fprintf('      %d comparaciones, 0 discrepancias -> las 3 filas son representativas\n', nChk);
assert(all(all(repAll == [1 2 3])), ...
    '6.10: los representantes derivados no son [1 2 3] en algun eje (%s).', mat2str(repAll));

% --- (1)-(4) los KPIs de cada eje contra las tablas PUBLICADAS de la memoria ------
%  Escritas aqui a mano a proposito: si se leyeran del .mat la comprobacion seria una
%  identidad y no comprobaria nada.
PUB = { [-11.61 -10.21  -6.48  -5.55  -4.97  -4.61  -4.51;
          -6.22  -4.56   0.82   2.57   3.80   4.56   4.74;
          -4.63  -2.79   3.82   6.04   7.25   7.74   7.86], ...
        [ -4.91  -4.80  -4.93  -5.21  -5.42  -5.76;
           3.93   4.29   3.97   3.25   2.80   2.04;
           7.56   7.83   7.41   6.58   5.58   4.56], ...
        [ -5.66  -7.41  -8.37  -9.91 -11.74 -13.19 -14.78 -14.93;
           2.29  -0.73  -2.17  -4.20  -6.38  -8.02  -9.74  -9.90;
           5.45   1.82   0.08  -2.39  -4.81  -6.56  -8.44  -8.56], ...
        [ -4.95  -4.95;   3.87   3.87;   7.23   7.23] };
NOM = {'(1) densidad orbital','(2) densidad de haces (var. B)', ...
       '(3) mascara de elevacion (T=66)','(4) operadores co-canal'};
for e = 1:4
    G = nan(3, numel(EJ(e).R));
    for r = 1:3
        G(r,:) = arrayfun(@(i) EJ(e).R{i}.K{repAll(e,r)}.viab.SINR_edge_p5, 1:numel(EJ(e).R));
    end
    d = max(max(abs(G - PUB{e})));
    chk(sprintf('6.10 %s SINR_edge_p5', NOM{e}), d, 0.006);
    assert(d <= 0.006, '6.10: el eje %s no reproduce la tabla publicada (%.4f dB).', NOM{e}, d);
    if e == 2
        % AVISO: no se esconde dentro del max. La ultima celda de reuse1 difiere mas
        % que el resto porque la TABLA PUBLICADA tenia un digito mal, no el .mat.
        % La referencia de arriba se deja EN EL VALOR PUBLICADO (-5.76) a proposito:
        % relajarla seria borrar la senal. El valor correcto segun el .mat es -5.75;
        % cuando la memoria se actualice, cambiar aqui PUB{2}(1,end) a -5.75 y esta
        % desviacion bajara a 0.0002 dB.
        dCell = abs(G(1,end) - PUB{2}(1,end));
        fprintf('      AVISO celda reuse1 @ dens=%.2f: publicado %.2f | .mat %.4f | dif %.4f dB\n', ...
            E310.dens(end), PUB{2}(1,end), G(1,end), dCell);
        fprintf('      -> pasa la tolerancia, pero es un DIGITO MAL EN LA TABLA PUBLICADA\n');
    end
end

% --- (5) la afirmacion que sostiene la figura, comprobada como tal ---------------
nCells = 0;  nViable = 0;
for e = 1:4
    for i = 1:numel(EJ(e).R)
        nCells  = nCells + 1;
        nViable = nViable + strcmp(EJ(e).R{i}.K{repAll(e,1)}.viab.verdict, 'viable');
    end
end
fprintf('  fila reuse1: %d celdas en los 4 ejes, de las cuales ''viable'' = %d\n', nCells, nViable);
chk('6.10 (5) reuse1 nunca viable', nViable, 0);
assert(nViable == 0, ...
    '6.10: reuse1 alcanza ''viable'' en %d celdas. Cae el argumento del capitulo.', nViable);
fprintf('  VERIFICACION OK (0..5) -> se promueve la figura\n');

% ------------------------------- dibujo -----------------------------------------
%  `lev` es el mapa que ya definen los bloques 6.5/6.7: se REUTILIZA en vez de
%  redeclararlo, para que los tres mapas del capitulo no puedan divergir.
assert(exist('lev','var') == 1 && isKey(lev,'viable') && isKey(lev,'inviable'), ...
    '6.10: no esta disponible el mapa `lev` de los bloques 6.5/6.7.');
CMAPV = [0.85 0.85 0.85; 0.80 0.20 0.20; 0.95 0.75 0.20; 0.20 0.65 0.30];
W22   = [100 100 760 560];              % 2x2: mas alto que W2/WL, que son de 1x2
ROWLB = {'reuse1','D = 3','D = 4'};

f = figure('Color','w','Visible','off','Position',W22);
for e = 1:4
    ax = subplot(2,2,e);
    nCol = numel(EJ(e).R);
    Mv = zeros(3, nCol);
    for r = 1:3
        for i = 1:nCol
            Mv(r,i) = lev(EJ(e).R{i}.K{repAll(e,r)}.viab.verdict);
        end
    end
    imagesc(Mv,[0 3]);
    colormap(ax, CMAPV);
    set(ax,'YTick',1:3,'YTickLabel',ROWLB);
    set(ax,'XTick',1:nCol,'XTickLabel',EJ(e).lab);
    if EJ(e).rot ~= 0, set(ax,'XTickLabelRotation',EJ(e).rot); end
    xlabel(EJ(e).xlab);
    % Abreviatura de 3 letras dentro de cada celda: sin ella el mapa solo se lee por
    % color y se pierde impreso en escala de grises (mismo criterio que 6.5/6.7).
    for r = 1:3
        for i = 1:nCol
            text(i, r, EJ(e).R{i}.K{repAll(e,r)}.viab.verdict(1:3), ...
                 'HorizontalAlignment','center', 'FontSize',11, ...
                 'Color','w', 'FontWeight','bold');
        end
    end
    panelOut(ax, EJ(e).tag);
end
bigfont(f); save_fig(f, STAGE, 'fig_6_10_mapa_veredictos_unificado'); close(f);
promote('fig_6_10_mapa_veredictos_unificado');

fprintf('\n[cap 6] terminado.\n');
