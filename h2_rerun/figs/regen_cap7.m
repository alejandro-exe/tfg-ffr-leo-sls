%% REGEN_CAP7   Las SIETE figuras del capitulo 7 (H2, FFR adaptativa), SIN TITULO.
%
%  RE-GRAFICADO PURO. Todo lo que se dibuja se LEE de los .mat: no se simula, no se
%  llama al motor (compute_sinr_ffr, compute_kpis, ffr_allocate, ffr_policy,
%  run_one_density, run_sweep_points) y no se modifica ningun .mat. La unica
%  aritmetica que se hace aqui es de PRESENTACION (cuartiles, medias por bloque) y
%  la reconstruccion de la regla INVERTIDA, que es la formula de ffr_policy aplicada
%  al alpha nominal ya almacenado -- no una reejecucion de la politica.
%
%  Estilo: plantilla de figuras/scripts/regen_cap6.m. Sin title() ni sgtitle (el
%  pie lo lleva el documento); ejes y leyendas en ASCII sin tildes; fuentes 11 pt;
%  300 dpi via save_fig; un bloque %% FIG 7.x por figura, de modo que cada figura
%  tenga su script y se puedan ANADIR bloques sin tocar los existentes.
%
%  VERIFICACION. Las referencias numericas estan ESCRITAS AQUI, no leidas del mismo
%  fichero que se grafica: comparar un .mat consigo mismo seria una identidad. Cada
%  figura comprueba las suyas ANTES de exportar y ABORTA si fallan. La unica
%  excepcion declarada es el percentil 95 del alpha oraculo, donde la referencia
%  externa usa OTRA convencion de percentil que prctile (ver FIG 7.04): se comprueba
%  contra las dos y se AVISA de la diferencia de convencion en vez de abortar.
%
%  FUENTES:
%    e5_adaptive_results.mat  E5 (T=1584): barrido de alpha, 4 casos, bins, alpha(t)
%    h2_anchors.mat           barrido de tau y anclas (contraste)
%    h2_oracle_v2.mat         frontera REFINADA y alpha ORACULO (v2 corregida)
%    h2_criterion.mat         tabla consolidada de Delta*SE_c/SE_b
%
%  NOMBRES HEREDADOS. Las figuras 7.01 y 7.05 se exportan tambien con el nombre que
%  tenian en run_e5_adaptive, pero con PREFIJO legacy_ y SIN numero de capitulo:
%  antes convivian en la misma carpeta dos ficheros distintos con el mismo prefijo
%  (fig_7_02_* y fig_7_03_*), y un \includegraphics mal autocompletado en Overleaf
%  coge la figura equivocada y COMPILA SIN ERROR. El prefijo legacy_ elimina la
%  colision de raiz. Lo mismo con el panel de cuartiles suelto de la 7.04.
%
%  SALIDA: h2_rerun/figs/*.png   LOG: h2_rerun/logs/regen_cap7.txt
%
%  Uso: run('.../h2_rerun/figs/regen_cap7.m')
%
%  TFG UC3M - Viabilidad de FFR adaptativa en constelaciones LEO.

clear; clc;

%% ---------------------------- rutas y log --------------------------------
HERE     = fileparts(mfilename('fullpath'));       % h2_rerun/figs
H2DIR    = fileparts(HERE);                        % h2_rerun
PROJROOT = fileparts(H2DIR);                       % raiz del proyecto
addpath(PROJROOT);                                 % save_fig vive en la raiz
FIGDIR = HERE;
LOGDIR = fullfile(H2DIR,'logs');
if ~exist(LOGDIR,'dir'), mkdir(LOGDIR); end
LOGF = fullfile(LOGDIR,'regen_cap7.txt');
diary off; if exist(LOGF,'file'), delete(LOGF); end; diary(LOGF);

bigfont = @(fh) set(findall(fh,'-property','FontSize'),'FontSize',11);
panel   = @(ax,t) text(ax, 0.02, 0.98, t, 'Units','normalized','FontSize',12, ...
            'FontWeight','bold','VerticalAlignment','top','HorizontalAlignment','left');
chk     = @(tag,d,tol) fprintf('  %-42s max|dif| = %.6f   (tol %.6f)\n', tag, d, tol);
W3 = [100 100 1080 340];                 % 3 paneles
WL = [100 100  920 460];                 % 1 panel con leyenda DENTRO del eje
W2 = [100 100 1120 420];                 % 2 paneles

fprintf('\n################ CAPITULO 7 (siete figuras) ################\n');
fprintf('salida: %s\n', FIGDIR);

%% ---------------------------- carga de fuentes ---------------------------
E5 = load(fullfile(PROJROOT,'e5_adaptive_results.mat'), ...
          'sw','GR','FR','alphas','KK','cases','alpha_star','DELTA', ...
          'elev_t','alp_t','ctxLite','bin','binName','ratio_bin','ratio_all');
OR = load(fullfile(H2DIR,'results','h2_oracle_v2.mat'), 'OB','TAUS','cfg');
CR = load(fullfile(H2DIR,'results','h2_criterion.mat'), 'H2');
AN = load(fullfile(H2DIR,'results','h2_anchors.mat'),   'tau_sweep','tau_list');

iR1 = 1; iST = 2; iAD = 3; iIN = 4;                % run_e5_adaptive.m:180
nm5 = cellfun(@(c) c.name, E5.cases, 'UniformOutput', false);
nA  = numel(E5.alphas);
% OB se guarda como CELL; se admite tambien struct array por si cambia el formato.
if iscell(OR.OB), O1 = OR.OB{1}; else, O1 = OR.OB(1); end
fprintf('esquemas E5: %s | Delta=%d | alpha* = %.2f\n', strjoin(nm5,' | '), E5.DELTA, E5.alpha_star);
fprintf('oraculo v2: tau = %s | rejilla refinada de %d alpha\n', ...
        mat2str(OR.TAUS), numel(O1.swR.alpha));

%% ========================================================================
%% FIG 7.01  alpha(t) y elevacion media del servidor (doble eje)
%%   Es la antigua fig_7_02_alpha_vs_elevacion, SOLO RENUMERADA. Se exporta con el
%%   nombre nuevo y TAMBIEN con el viejo, para no romper referencias existentes.
%% ========================================================================
fprintf('\n--- FIG 7.01  alpha(t) vs elevacion ---\n');
th = (E5.ctxLite.tvec - E5.ctxLite.tvec(1))/60;        % min
rA = [min(E5.alp_t)  max(E5.alp_t)];
rE = [min(E5.elev_t) max(E5.elev_t)];
cc = corr(E5.alp_t(:), E5.elev_t(:));
fprintf('  alpha(t) en [%.4f, %.4f] (ref [0.4676, 0.6725])\n', rA);
fprintf('  elev(t)  en [%.2f, %.2f] deg (ref [52.23, 85.53])\n', rE);
fprintf('  corr(alpha, elev) = %+.6f (ref +1.000000)\n', cc);
d1a = max(abs(rA - [0.4676 0.6725]));  chk('7.01 recorrido de alpha', d1a, 0.0002);
d1b = max(abs(rE - [52.23 85.53]));    chk('7.01 recorrido de elevacion', d1b, 0.01);
d1c = abs(cc - 1);                     chk('7.01 corr(alpha,elev)', d1c, 5e-4);
assert(d1a <= 0.0002 && d1b <= 0.01 && d1c <= 5e-4, '7.01: E5 no reproduce las referencias.');
assert(abs(E5.alpha_star - 0.55) <= 1e-12, '7.01: alpha* no es 0.55.');

f = figure('Color','w','Visible','off');
yyaxis left
plot(th, E5.elev_t, 'LineWidth',1.6);
ylabel('elevacion media del servidor [deg]');
yyaxis right
plot(th, E5.alp_t, 'LineWidth',1.8); hold on;
plot(th, E5.alpha_star*ones(size(th)), ':', 'LineWidth',1.6);
ylabel('\alpha, fraccion de banda del INTERIOR [-]'); ylim([0 1]);
grid on; box on; xlabel('tiempo [min]');
legend({'elevacion media del servidor','\alpha adaptativa', ...
        sprintf('\\alpha estatica = %.2f', E5.alpha_star)}, 'Location','best');
bigfont(f);
save_fig(f, FIGDIR, 'fig_7_01_alpha_vs_elevacion');
% Nombre heredado, con prefijo legacy_ y SIN numero de capitulo: como
% fig_7_02_alpha_vs_elevacion colisionaba con fig_7_02_pareto_refinado, un
% \includegraphics autocompletado en Overleaf podia coger la otra y compilar bien.
save_fig(f, FIGDIR, 'legacy_e5_c_alpha_vs_elevacion');
close(f);
fprintf('  exportada como fig_7_01_alpha_vs_elevacion (+ heredada legacy_e5_c_...)\n');

%% ========================================================================
%% FIG 7.02  Plano (R_p5 global, R_p5 borde) con la frontera REFINADA
%%   Barrido grueso (19 puntos) + frontera refinada del oraculo v2 (paso 0.01 en la
%%   ventana de alpha*) + los cuatro puntos de E5.
%%
%%   LA GANANCIA NO SE DIBUJA. Son -0.094 Mbps sobre un eje que llega a 16: el trazo
%%   es INVISIBLE, y una leyenda que anuncia una linea que no se ve es peor que no
%%   ponerla. La cifra va en el PIE de figura. Se sigue verificando aqui (y se
%%   imprime en el log) para que quien escriba el pie tenga el numero exacto: es la
%%   medida contra la frontera REFINADA (-0.094), no contra la cuerda gruesa
%%   (-0.063), que subestima la frontera y por tanto INFLA la ganancia.
%% ========================================================================
fprintf('\n--- FIG 7.02  Pareto refinado ---\n');
gG = cellfun(@(K) K.all.R_p5_Mbps,  E5.KK);
gE = cellfun(@(K) K.edge.R_p5_Mbps, E5.KK);
ref_pts = [5.948890 5.589765;      % reuso-1
           7.446712 7.430042;      % estatica alpha* = 0.55
           7.184932 8.267286;      % adaptativa
           7.222372 8.066918];     % invertida
got_pts = [gG(iR1) gE(iR1); gG(iST) gE(iST); gG(iAD) gE(iAD); gG(iIN) gE(iIN)];
for r = 1:4
    fprintf('  %-16s (%.6f, %.6f)  ref (%.6f, %.6f)\n', nm5{r}, ...
        got_pts(r,1), got_pts(r,2), ref_pts(r,1), ref_pts(r,2));
end
d2a = max(abs(got_pts(:) - ref_pts(:)));  chk('7.02 los cuatro puntos de E5', d2a, 1e-6);
fprintf('  frontera GRUESA: %d no dominados de %d (ref 11 de 19)\n', numel(E5.GR.fx), nA);
d2b = abs(E5.GR.fx(1) - 0.734820);   chk('7.02 GR.fx primero', d2b, 1e-6);
d2c = abs(E5.GR.fx(end) - 7.446712); chk('7.02 GR.fx ultimo',  d2c, 1e-6);
gref = O1.gain_nominal_ref;
fprintf('  ganancia contra la frontera REFINADA = %+.6f Mbps -> PARA EL PIE: %+.3f\n', gref, gref);
fprintf('  ganancia contra la cuerda GRUESA     = %+.6f Mbps (no se cita)\n', O1.gain_nominal);
d2d = abs(round(gref,3) - (-0.094));  chk('7.02 ganancia refinada = -0.094', d2d, 1e-9);
% La leyenda va DENTRO, en 'northeast': hay que comprobar que ese cuadrante esta
% realmente vacio y no se estan tapando datos. Caja generosa (x >= 60% del rango,
% y >= 70%), que es mas de lo que ocupa la leyenda a 9 pt.
xall = [E5.sw.R_p5(:); O1.GRstatR.fx(:); got_pts(:,1)];
yall = [E5.sw.R_p5_edge(:); O1.GRstatR.fy(:); got_pts(:,2)];
xlim2 = min(xall) + 0.60*(max(xall)-min(xall));
ylim2 = 0.70*18;                                   % el ylim que se fija abajo
nBox  = sum(xall >= xlim2 & yall >= ylim2);
fprintf('  puntos en la caja de la leyenda (x>=%.2f, y>=%.2f): %d (ref 0)\n', xlim2, ylim2, nBox);
assert(d2a <= 1e-6 && numel(E5.GR.fx) == 11 && nA == 19 && ...
       d2b <= 1e-6 && d2c <= 1e-6 && d2d <= 1e-9, '7.02: no reproduce las referencias.');
assert(nBox == 0, '7.02: el cuadrante superior derecho ya NO esta vacio; la leyenda taparia datos.');

f = figure('Color','w','Visible','off','Position',WL);
plot(E5.sw.R_p5, E5.sw.R_p5_edge, 'o', 'Color',[.55 .55 .55], 'LineWidth',1.2, ...
     'MarkerSize',6, 'MarkerFaceColor','w', ...
     'DisplayName','FFR estatica, barrido grueso (19 \alpha)'); hold on;
plot(O1.GRstatR.fx, O1.GRstatR.fy, 'k-', 'LineWidth',2.2, ...
     'DisplayName','frontera de Pareto REFINADA (paso 0.01)');
for ia = 1:2:nA
    text(E5.sw.R_p5(ia), E5.sw.R_p5_edge(ia), sprintf('  %.2f', E5.alphas(ia)), ...
         'FontSize',9, 'Color',[.45 .45 .45]);
end
plot(gG(iR1), gE(iR1), 'ks', 'MarkerSize',10, 'MarkerFaceColor',[.8 .8 .8], ...
     'DisplayName','reuso-1');
plot(gG(iST), gE(iST), 'bd', 'MarkerSize',10, 'MarkerFaceColor','b', ...
     'DisplayName',sprintf('FFR estatica (\\alpha*=%.2f)', E5.alpha_star));
plot(gG(iAD), gE(iAD), 'rp', 'MarkerSize',16, 'MarkerFaceColor','r', ...
     'DisplayName','FFR adaptativa (H2)');
plot(gG(iIN), gE(iIN), 'v', 'Color',[.85 .5 .1], 'MarkerSize',9, ...
     'MarkerFaceColor',[.85 .5 .1], 'DisplayName','adaptativa INVERTIDA');
grid on; box on;
xlabel('R_{p5} global [Mbps]');
ylabel('R_{p5} de los usuarios de BORDE [Mbps]');
% LEYENDA DENTRO DEL EJE, en 'northeast'. Los datos ocupan la ANTI-diagonal, luego
% el cuadrante superior derecho esta vacio y ahi no tapa nada; se comprueba abajo
% que no hay ningun punto dentro de la caja de la leyenda. Fuera del eje NO vale:
% a \linewidth en la memoria la leyenda se comeria ~40% del ancho del grafico.
ylim([0 18]);                                          % hueco para la leyenda
lg = legend('Location','northeast');
bigfont(f);
set(lg,'FontSize',9);                                  % DESPUES de bigfont (que pone 11)
save_fig(f, FIGDIR, 'fig_7_02_pareto_refinado'); close(f);
fprintf('  exportada como fig_7_02_pareto_refinado (leyenda DENTRO, sin linea de ganancia)\n');

%% ========================================================================
%% FIG 7.03  Criterio consolidado  rho = Delta*SE_c/SE_b, por TRES vias
%%   (a) E5 por tramo de elevacion, (b) barrido de densidad frente a T con Delta=3
%%   y Delta=4, (c) OneWeb HEXAGONAL. La reticula 4x4 se EXCLUYE a proposito: es
%%   material del capitulo 8 (topologia), no de la refutacion de H2.
%% ========================================================================
fprintf('\n--- FIG 7.03  criterio consolidado ---\n');
H = CR.H2;
% (a) E5
selA  = H.fuente == "E5";
casA  = H.caso(selA);   rhoA = H.crit(selA);
ordA  = [find(contains(casA,'25-35')), find(contains(casA,'35-60')), ...
         find(contains(casA,'60-90')), find(contains(casA,'global'))];
rhoA  = rhoA(ordA);
refA  = [NaN 0.777904 0.862729 0.855612];
fprintf('  (a) E5 por tramo + global: %s\n', mat2str(round(rhoA(:).',6)));
d3a = max(abs(rhoA(2:4) - refA(2:4).'));
chk('7.03a rho de E5', d3a, 1e-6);
assert(isnan(rhoA(1)) && d3a <= 1e-6, '7.03a: rho de E5 no coincide (o el tramo vacio no es NaN).');

% (b) barrido de densidad (solo FFR estatica; la adaptativa da el MISMO cociente)
selD = H.fuente == "densidad";
Td   = H.T(selD);  Dd = H.Delta(selD);  rd = H.crit(selD);  cd = H.caso(selD);
isSt = contains(cd,'FFR (D=');                    % excluye FFR-ADAPTATIVA
T3 = Td(isSt & Dd==3);  r3 = rd(isSt & Dd==3);
T4 = Td(isSt & Dd==4);  r4 = rd(isSt & Dd==4);
[T3,i3] = sort(T3); r3 = r3(i3);
[T4,i4] = sort(T4); r4 = r4(i4);
ref3 = [0.885 0.874 0.853 0.837 0.856 0.877 0.891];
ref4 = [0.846 0.832 0.843 0.883 0.949 0.997 1.020];
refT = [66 300 600 1000 1584 2500 4000];
fprintf('  (b) Delta=3: %s\n', mat2str(round(r3(:).',6)));
fprintf('  (b) Delta=4: %s\n', mat2str(round(r4(:).',6)));
% TOLERANCIA 1e-3, y hay que justificarla: las referencias de este bloque vienen
% dadas a TRES decimales, luego su propio redondeo vale ya 5e-4 y no tiene sentido
% exigir mas. Se comparan los valores CRUDOS contra la referencia (no los valores
% redondeados contra ella, que es lo que introduce un segundo redondeo). El punto
% mas ajustado es T=300 con Delta=4: el valor real es 0.831499 y la referencia dice
% 0.832, es decir 0.831499 -> 0.8315 -> 0.832, redondeado DOS veces; la desviacion
% resultante (5.0e-4) es exactamente el ancho del redondeo, no un desajuste de datos.
% Estos mismos valores ya estan verificados a 1e-9 contra el .mat en h2_criterion.
TOL3 = 1e-3;
d3b1 = max(abs(r3(:).' - ref3));  chk('7.03b rho Delta=3 vs T (ref a 3 dec)', d3b1, TOL3);
d3b2 = max(abs(r4(:).' - ref4));  chk('7.03b rho Delta=4 vs T (ref a 3 dec)', d3b2, TOL3);
d3b3 = max(abs(T3(:).' - refT));
[wm, iwm] = max(abs(r4(:).' - ref4));
fprintf('  punto mas ajustado de Delta=4: T=%d, valor %.6f vs referencia %.3f (dif %.1e)\n', ...
    T4(iwm), r4(iwm), ref4(iwm), wm);
assert(d3b1 <= TOL3 && d3b2 <= TOL3 && d3b3 == 0 && isequal(T3,T4), ...
    '7.03b: el barrido de densidad no coincide.');

% (c) OneWeb HEXAGONAL (la reticula 4x4 se excluye: capitulo 8)
selW = H.arquitectura == "OneWeb-hex" & contains(H.caso,'ffr(') & ~contains(H.caso,'adapt');
Dw   = H.Delta(selW);  rw = H.crit(selW);
[Dw,iw] = sort(Dw);    rw = rw(iw);
refW = [0.800421 0.929638];
fprintf('  (c) OneWeb hex Delta=%s: %s\n', mat2str(Dw(:).'), mat2str(round(rw(:).',6)));
d3c = max(abs(rw(:).' - refW));  chk('7.03c rho OneWeb hex', d3c, 1e-6);
assert(d3c <= 1e-6 && isequal(Dw(:).',[3 4]), '7.03c: OneWeb hexagonal no coincide.');
nRect = sum(H.arquitectura == "OneWeb-rect");
fprintf('  (c) filas de la reticula 4x4 EXCLUIDAS a proposito: %d\n', nRect);

f = figure('Color','w','Visible','off','Position',W3);
% (a)
ax = subplot(1,3,1);
bar(1:4, [0 rhoA(2) rhoA(3) rhoA(4)], 0.6, 'FaceColor',[.30 .50 .75]); hold on;
yline(1.0,'k--','LineWidth',1.6);
text(1, 0.05, 'vacio', 'Rotation',90, 'FontSize',10, 'Color',[.4 .4 .4], ...
     'HorizontalAlignment','left');
set(gca,'XTick',1:4,'XTickLabel',{'25-35','35-60','60-90','global'});
ylim([0 1.35]); grid on; box on;
xlabel('tramo de elevacion [deg]'); ylabel('\rho = \Delta\cdotSE_c/SE_b [-]');
panel(ax,'(a)');
% (b)
ax = subplot(1,3,2);
semilogx(T3, r3, 'o-', 'LineWidth',1.8, 'MarkerFaceColor','w', ...
         'DisplayName','\Delta = 3'); hold on;
semilogx(T4, r4, 's-', 'LineWidth',1.8, 'MarkerFaceColor','w', ...
         'DisplayName','\Delta = 4');
yline(1.0,'k--','LineWidth',1.6,'HandleVisibility','off');
ylim([0 1.35]); grid on; box on;
xlabel('satelites de la constelacion T [-]'); ylabel('\rho [-]');
legend('Location','southeast');
panel(ax,'(b)');
% (c)
ax = subplot(1,3,3);
bar(1:2, rw(:).', 0.5, 'FaceColor',[.85 .50 .15]); hold on;
yline(1.0,'k--','LineWidth',1.6);
set(gca,'XTick',1:2,'XTickLabel',{'\Delta = 3','\Delta = 4'});
ylim([0 1.35]); grid on; box on;
xlabel('OneWeb hexagonal (19 haces)'); ylabel('\rho [-]');
panel(ax,'(c)');
bigfont(f); save_fig(f, FIGDIR, 'fig_7_03_criterio_consolidado'); close(f);
fprintf('  exportada como fig_7_03_criterio_consolidado\n');

%% ========================================================================
%% FIG 7.04  alpha de las tres politicas: (a) por CUARTIL de elevacion del
%%           servidor y (b) DISTRIBUCION sobre los 121 instantes (tau = 50%)
%%
%%   Tres series: regla nominal, regla INVERTIDA y alpha ORACULO. La invertida NO
%%   esta almacenada, asi que se reconstruye con la formula de ffr_policy aplicada
%%   al alpha nominal YA GUARDADO (aritmetica de presentacion, no reejecucion):
%%       w = (a_nom - a_min)/(a_max - a_min)  ->  a_inv = a_max + (a_min-a_max)*w
%%   Como en esta configuracion alpha_min + alpha_max = 1, equivale a 1 - a_nom, y
%%   asi se comprueba (control a 1e-16).
%%
%%   (a) Cuartiles: instantes ORDENADOS por elevacion, en bloques de igual numero
%%       de instantes (30/30/30/31), con el resto en el ULTIMO. Cada entrada de la
%%       leyenda lleva su ganancia sobre la frontera REFINADA, que es lo que
%%       convierte la figura en un argumento y no en tres barras: el alpha que
%%       DECRECE con la elevacion es el unico que gana.
%%   (b) Distribucion de alpha(t). Enseña dos cosas que los cuartiles promedian y
%%       por tanto esconden: (i) la masa del oraculo esta ENTERA por debajo de la
%%       de las dos reglas, y (ii) el oraculo baja tres instantes hasta el SUELO
%%       alpha = 0.05, un valor que ninguna regla de elevacion puede alcanzar
%%       porque su recorrido esta acotado por las anclas.
%%
%%   La ganancia de la INVERTIDA no esta guardada en ningun .mat (el oraculo solo
%%   almacena la suya y la nominal), asi que se obtiene interpolando la frontera
%%   REFINADA ya almacenada (O1.GRstatR) en su R_p5 global. Es el mismo calculo que
%%   hace frontier_gain, y se VALIDA reproduciendo con el las dos ganancias que si
%%   estan guardadas antes de aplicarlo a la invertida.
%% ========================================================================
fprintf('\n--- FIG 7.04  alpha: cuartiles y distribucion ---\n');
elO = O1.elevMean(:).';
aNom = O1.alpha_nominal(:).';
aOra = O1.alpha_oracle(:).';
Ntq  = numel(elO);
aMin = OR.cfg.ffr.alpha_min;  aMax = OR.cfg.ffr.alpha_max;
w    = (aNom - aMin) / (aMax - aMin);
aInv = aMax + (aMin - aMax) * w;
fprintf('  anclas de alpha: [%.2f, %.2f] | w en [%.4f, %.4f]\n', aMin, aMax, min(w), max(w));
assert(all(w >= -1e-12 & w <= 1+1e-12), '7.04: w fuera de [0,1]; la reconstruccion no vale.');
% Comprobacion cruzada: con a_min + a_max = 1 la invertida es exactamente 1 - a_nom.
if abs((aMin + aMax) - 1) <= 1e-12
    fprintf('  control: max|a_inv - (1 - a_nom)| = %.3e\n', max(abs(aInv - (1 - aNom))));
end

qsz = repmat(floor(Ntq/4),1,4);  qsz(4) = Ntq - sum(qsz(1:3));   % 30/30/30/31
[~, ord] = sort(elO, 'ascend');
Q = struct('nom',nan(1,4),'inv',nan(1,4),'ora',nan(1,4),'elev',nan(1,4),'n',qsz);
p = 0;
for q = 1:4
    idx = ord(p+1 : p+qsz(q));  p = p + qsz(q);
    Q.nom(q)  = mean(aNom(idx));
    Q.inv(q)  = mean(aInv(idx));
    Q.ora(q)  = mean(aOra(idx));
    Q.elev(q) = mean(elO(idx));
end
% Version de la referencia EXTERNA: bloques de 30 descartando el instante 121.
ord30 = ord(1:4*floor(Ntq/4));
oraRef30 = arrayfun(@(q) mean(aOra(ord30((q-1)*30+1 : q*30))), 1:4);

fprintf('  tamanos de cuartil: %s (%d instantes)\n', mat2str(qsz), Ntq);
fprintf('  elevacion media por cuartil : %s deg\n', mat2str(round(Q.elev,2)));
fprintf('  alpha NOMINAL   por cuartil : %s\n', mat2str(round(Q.nom,4)));
fprintf('  alpha INVERTIDA por cuartil : %s\n', mat2str(round(Q.inv,4)));
fprintf('  alpha ORACULO   por cuartil : %s   (30/30/30/31)\n', mat2str(round(Q.ora,4)));
fprintf('  alpha ORACULO   por cuartil : %s   (30x4, descartando el ultimo)\n', ...
        mat2str(round(oraRef30,4)));

% Monotonias exigidas (bloqueantes)
cN = corr(Q.nom(:), Q.elev(:));  cI = corr(Q.inv(:), Q.elev(:));
fprintf('  corr(alpha_nominal, elev) por cuartil   = %+.6f (ref +1, creciente)\n', cN);
fprintf('  corr(alpha_invertida, elev) por cuartil = %+.6f (ref -1, decreciente)\n', cI);
assert(all(diff(Q.nom) > 0) && abs(cN - 1) <= 1e-6, ...
    '7.04: alpha nominal no CRECE con la elevacion.');
assert(all(diff(Q.inv) < 0) && abs(cI + 1) <= 1e-6, ...
    '7.04: alpha invertida no DECRECE con la elevacion.');

% Cuartiles del oraculo: referencia AHORA BLOQUEANTE a 1e-4 (con la v1 rota habia
% que avisar en vez de abortar; la v2 corregida reproduce la referencia externa).
refOra = [0.4439 0.4328 0.3890 0.3619];
d4a = max(abs(Q.ora - refOra));
d4b = max(abs(oraRef30 - refOra));
chk('7.04 alpha oraculo por cuartil (30/30/30/31)', d4a, 1e-4);
assert(d4a <= 1e-4, '7.04: los cuartiles del alpha oraculo no reproducen la referencia.');
% El reparto 30x4 (descartando el instante 121) se imprime solo como CONTRASTE: no
% es el que se dibuja ni el de la referencia, y su Q4 se desvia 4.3e-4 por prescindir
% de una muestra. Sirve para ver que el resultado no depende del reparto elegido.
fprintf('  contraste, reparto 30x4 (no es el dibujado): max|dif| vs ref = %.6f\n', d4b);

%% ---- referencias de la DISTRIBUCION de alpha (panel b) ------------------
mNom = mean(aNom);  mInv = mean(aInv);  mOra = mean(aOra);
rNom = [min(aNom) max(aNom)];  rInv = [min(aInv) max(aInv)];  rOra = [min(aOra) max(aOra)];
nFloor = sum(aOra <= 0.06);
fprintf('  alpha medio    : nominal %.6f | invertida %.6f | oraculo %.6f\n', mNom, mInv, mOra);
fprintf('  recorrido nom  : [%.4f, %.4f]   (ref [0.4676, 0.6725])\n', rNom);
fprintf('  recorrido inv  : [%.4f, %.4f]   (ref [0.3275, 0.5324])\n', rInv);
fprintf('  recorrido ora  : [%.4f, %.4f]   (ref [0.0500, 0.5096])\n', rOra);
fprintf('  instantes del oraculo con alpha <= 0.06: %d de %d (ref 3 de 121)\n', nFloor, Ntq);
d4c = max(abs([mNom mInv mOra] - [0.568496 0.431504 0.406544]));
d4d = max(abs([rNom rInv rOra] - [0.4676 0.6725 0.3275 0.5324 0.0500 0.5096]));
chk('7.04 alpha medio de las tres politicas', d4c, 1e-4);
chk('7.04 recorridos de alpha',               d4d, 1e-4);
assert(d4c <= 1e-4 && d4d <= 1e-4, '7.04: medias o recorridos de alpha no coinciden.');
assert(nFloor == 3 && Ntq == 121, '7.04: no hay 3 instantes en el suelo, o Nt no es 121.');

% PERCENTILES DEL ALPHA ORACULO. La referencia externa NO usa la convencion de
% prctile: en p95 da 0.4961 y prctile da 0.496368. No es un desajuste de datos sino
% de DEFINICION de percentil -- con n=121, 0.95*(n-1) = 114 cae EXACTAMENTE sobre
% una muestra ordenada (la 115a, 0.496140), mientras que prctile usa (i-0.5)/n e
% interpola entre la 115a y la 116a. Se comprueban las DOS y se aborta solo si
% falla la convencion que la referencia usa; la de prctile (que es la del motor,
% compute_kpis) se imprime con su desviacion para que quede constancia.
qlev  = [5 25 50 75 95];
refPct = [0.3500 0.3750 0.4050 0.4500 0.4961];
pctML  = prctile(aOra, qlev);
aS = sort(aOra);  pos = qlev/100*(Ntq-1);  loI = floor(pos);  frI = pos - loI;
pctLIN = aS(loI+1).*(1-frI) + aS(min(loI+2,Ntq)).*frI;
fprintf('  percentiles %s del alpha oraculo:\n', mat2str(qlev));
fprintf('    interpolacion lineal (convencion de la referencia): %s\n', mat2str(round(pctLIN,6)));
fprintf('    prctile de MATLAB   (convencion de compute_kpis)  : %s\n', mat2str(round(pctML,6)));
d4e = max(abs(pctLIN - refPct));
d4f = max(abs(pctML  - refPct));
chk('7.04 percentiles (convencion de la referencia)', d4e, 1e-4);
chk('7.04 percentiles (prctile, solo informativo)',   d4f, 1e-4);
if d4f > 1e-4
    fprintf(['  AVISO (no aborta): con prctile el p95 sale %.6f frente a %.4f de la\n' ...
             '         referencia (dif %.1e). Es diferencia de CONVENCION de percentil,\n' ...
             '         no de datos: los otros cuatro coinciden con las dos convenciones.\n'], ...
        pctML(end), refPct(end), d4f);
end
assert(d4e <= 1e-4, '7.04: los percentiles del alpha oraculo no reproducen la referencia.');

% Punto estatico de contraste del panel (b): el alpha CONSTANTE mas proximo a la
% media del oraculo que esta medido en el barrido, y el objetivo que el oraculo
% respeta. Es lo que enseña que el oraculo no es "bajar alpha": con alpha bajo pero
% CONSTANTE el percentil global se hunde, porque los usuarios de CENTRO se quedan
% sin banda en TODOS los instantes.
ia40 = find(abs(E5.alphas - 0.40) <= 1e-12, 1);
assert(~isempty(ia40), '7.04: alpha = 0.40 no esta en el barrido de E5.');
R40  = E5.sw.R_p5(ia40);
targ = O1.target;
fprintf('  estatica alpha = 0.40 -> R_p5 global = %.6f (ref 5.878557)\n', R40);
fprintf('  objetivo del oraculo (R_p5 global de la adaptativa) = %.6f (ref 7.184932)\n', targ);
d4g = abs(R40 - 5.878557);  chk('7.04 estatica alpha=0.40', d4g, 1e-4);
d4h = abs(targ - 7.184932);  chk('7.04 objetivo del oraculo', d4h, 1e-4);
assert(d4g <= 1e-4 && d4h <= 1e-4, '7.04: el contraste estatico/objetivo no coincide.');

%% ---- ganancias sobre la frontera REFINADA (leyenda del panel a) ---------
% La de la invertida NO esta guardada: se interpola la frontera refinada, con el
% metodo VALIDADO contra las dos ganancias que si estan guardadas.
[fxs, ifx] = sort(O1.GRstatR.fx(:));  fys = O1.GRstatR.fy(:);  fys = fys(ifx);
gain_ref = @(x,y) y - interp1(fxs, fys, x, 'linear');
vN = gain_ref(gG(iAD), gE(iAD));                     % debe dar gain_nominal_ref
vO = gain_ref(O1.p5g_oracle, O1.p5e_oracle);         % debe dar gain_oracle_ref
dV = max(abs([vN vO] - [O1.gain_nominal_ref O1.gain_oracle_ref]));
chk('7.04 interpolacion de la frontera reproduce las ganancias guardadas', dV, 1e-6);
assert(dV <= 1e-6, '7.04: la interpolacion de la frontera NO reproduce lo guardado.');
gNom = O1.gain_nominal_ref;
gOra = O1.gain_oracle_ref;
gInv = gain_ref(gG(iIN), gE(iIN));                   % la que no esta guardada
fprintf('  ganancia sobre la frontera REFINADA: nominal %+.6f | invertida %+.6f | oraculo %+.6f\n', ...
        gNom, gInv, gOra);
d4i = max(abs(round([gNom gInv gOra],3) - [-0.094 -0.224 0.456]));
chk('7.04 ganancias redondeadas a 3 decimales', d4i, 1e-9);
assert(d4i <= 1e-9, '7.04: las ganancias no son -0.094 / -0.224 / +0.456.');

%% ---- dibujo ------------------------------------------------------------
% Etiquetas de tick en UNA sola linea: un '\n' dentro de un XTickLabel de cell
% array NO crea dos renglones, MATLAB lo parte en dos ticks distintos y la figura
% sale con los rangos desalineados y sin Q3/Q4.
qlab = cell(1,4);
for q = 1:4
    idxq = ord(sum(qsz(1:q-1))+1 : sum(qsz(1:q)));
    qlab{q} = sprintf('Q%d (%.0f-%.0f)', q, min(elO(idxq)), max(elO(idxq)));
end
C1 = [0.00 0.447 0.741];      % nominal
C2 = [0.85 0.325 0.098];      % invertida
C3 = [0.929 0.694 0.125];     % oraculo
legGain = {sprintf('regla nominal: %+.3f Mbps', gNom), ...
           sprintf('regla INVERTIDA: %+.3f Mbps', gInv), ...
           sprintf('\\alpha ORACULO (sin regla): %+.3f Mbps', gOra), ...
           sprintf('\\alpha* estatica = %.2f', E5.alpha_star)};

f = figure('Color','w','Visible','off','Position',W2);

% ---- panel (a): alpha medio por cuartil de elevacion --------------------
ax = subplot(1,2,1);
draw_quartiles(Q, qlab, E5.alpha_star, [C1;C2;C3]);
ylim([0 1.18]);
lg = legend(legGain, 'Location','northeast');
lg.Title.String = 'ganancia sobre la frontera refinada';
panel(ax,'(a)');

% ---- panel (b): distribucion de alpha sobre los 121 instantes ----------
ax = subplot(1,2,2);
edges = 0:0.025:0.7;                                  % <-- DECISION (bin de 0.025)
histogram(aNom, edges, 'FaceColor',C1, 'FaceAlpha',0.60, 'EdgeColor',[.25 .25 .25], ...
          'DisplayName', sprintf('regla nominal (media %.3f)', mNom)); hold on;
histogram(aInv, edges, 'FaceColor',C2, 'FaceAlpha',0.60, 'EdgeColor',[.25 .25 .25], ...
          'DisplayName', sprintf('regla INVERTIDA (media %.3f)', mInv));
histogram(aOra, edges, 'FaceColor',C3, 'FaceAlpha',0.60, 'EdgeColor',[.25 .25 .25], ...
          'DisplayName', sprintf('\\alpha ORACULO (media %.3f)', mOra));
ymax = 0;
for h = findobj(ax,'Type','Histogram').'
    ymax = max(ymax, max(h.Values));
end
ytop = ceil(ymax*1.55/5)*5;                            % hueco para leyenda y notas
yLin = ytop*0.78;                                      % las verticales paran BAJO la leyenda
% Medias de cada politica (sin entrada de leyenda: ya van en el texto de cada una)
for k = 1:3
    mk = [mNom mInv mOra]; ck = [C1;C2;C3];
    plot([mk(k) mk(k)], [0 yLin], '--', 'Color',ck(k,:), 'LineWidth',1.8, ...
         'HandleVisibility','off');
end
plot([E5.alpha_star E5.alpha_star], [0 yLin], 'k:', 'LineWidth',1.8, ...
     'DisplayName', sprintf('\\alpha* estatica = %.2f', E5.alpha_star));
% El suelo: 3 instantes en alpha = 0.05. Es una barra corta y aislada, asi que se
% señala explicitamente; es justo lo que ninguna regla acotada por anclas alcanza.
% La banda 0.075 < alpha < 0.33 esta practicamente vacia (el p5 del oraculo es
% 0.35), luego el rotulo cabe ahi sin tapar ninguna barra.
plot([0.0625 0.072], [nFloor+0.4 nFloor+1.4], 'k-', 'LineWidth',1.0, 'HandleVisibility','off');
text(0.075, nFloor+1.5, sprintf('%d instantes en\nel suelo \\alpha = %.2f', nFloor, min(aOra)), ...
     'FontSize',9, 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
% Rotulo del contraste estatico. Va en la MITAD IZQUIERDA y a media altura: arriba
% lo tapa la leyenda y abajo lo tapa el rotulo del suelo. La banda alpha < 0.32 esta
% vacia en toda la altura, luego no oculta ninguna barra.
text(0.02, 0.55, sprintf(['media del \\alpha oraculo = %.2f,\n' ...
     'pero un \\alpha CONSTANTE = 0.40\n' ...
     'da R_{p5} global = %.2f, frente\n' ...
     'al objetivo %.2f'], mOra, R40, targ), 'Units','normalized', ...
     'FontSize',9, 'VerticalAlignment','top', 'HorizontalAlignment','left', ...
     'BackgroundColor','w', 'EdgeColor',[.6 .6 .6], 'Margin',3);
xlim([0 0.7]); ylim([0 ytop]); grid on; box on;
xlabel('\alpha, fraccion de banda del INTERIOR [-]');
ylabel('instantes [-]');
lg2 = legend('Location','northeast');
panel(ax,'(b)');

bigfont(f);
set([lg lg2],'FontSize',8.5);                          % DESPUES de bigfont
save_fig(f, FIGDIR, 'fig_7_04_alpha_oraculo'); close(f);
fprintf('  exportada como fig_7_04_alpha_oraculo (dos paneles)\n');

% Panel (a) suelto, con la leyenda ORIGINAL (sin ganancias): es el fichero que
% habia antes de partir la figura en dos. Se conserva con prefijo legacy_ para no
% colisionar con fig_7_04_alpha_oraculo, y regenerable en vez de huerfano.
f = figure('Color','w','Visible','off');
draw_quartiles(Q, qlab, E5.alpha_star, [C1;C2;C3]);
ylim([0 0.95]);
legend({'regla nominal (\alpha crece con la elevacion)', ...
        'regla INVERTIDA', ...
        '\alpha ORACULO (sin regla)', ...
        sprintf('\\alpha* estatica = %.2f', E5.alpha_star)}, 'Location','northeast');
bigfont(f); save_fig(f, FIGDIR, 'legacy_fig_7_04_cuartiles'); close(f);
fprintf('  exportada la heredada legacy_fig_7_04_cuartiles (solo el panel de cuartiles)\n');

%% ========================================================================
%% FIG 7.05  R_p5 de BORDE por tramo de elevacion
%%   Es la antigua fig_7_03_rp5_por_tramo_elevacion, SOLO RENUMERADA. Se exporta
%%   con el nombre nuevo y TAMBIEN con el viejo.
%% ========================================================================
fprintf('\n--- FIG 7.05  R_p5 de borde por tramo ---\n');
B = E5.bin.R_p5_edge;
fprintf('  tramos: %s\n', strjoin(E5.binName,' | '));
for c = [iR1 iST iAD]
    fprintf('  %-16s %s\n', nm5{c}, mat2str(round(B(c,:),4)));
end
fprintf('  tramo "%s" VACIO: %d (ref 1)\n', E5.binName{1}, all(isnan(B(:,1))));
g2 = 100*(B(iAD,2)/B(iST,2) - 1);
g3 = 100*(B(iAD,3)/B(iST,3) - 1);
fprintf('  adaptativa vs estatica: %+.2f %% (35-60) | %+.2f %% (60-90)  ref +15.6 / -5.6\n', g2, g3);
d5a = max(abs([g2 g3] - [15.6 -5.6]));  chk('7.05 ganancias por tramo', d5a, 0.06);
assert(all(isnan(B(:,1))) && d5a <= 0.06, '7.05: no reproduce las referencias.');

f = figure('Color','w','Visible','off');
bar(B([iR1 iST iAD],:).');
set(gca,'XTickLabel',E5.binName);
grid on; box on;
ylabel('R_{p5} de BORDE [Mbps]');
xlabel('tramo de elevacion del servidor [deg]');
legend(nm5([iR1 iST iAD]),'Location','northwest');
bigfont(f);
save_fig(f, FIGDIR, 'fig_7_05_rp5_por_tramo');
% Nombre heredado con prefijo legacy_ y SIN numero: fig_7_03_rp5_por_tramo_elevacion
% colisionaba con fig_7_03_criterio_consolidado (misma trampa que en la 7.01).
save_fig(f, FIGDIR, 'legacy_e5_d_rp5_por_tramo');
close(f);
fprintf('  exportada como fig_7_05_rp5_por_tramo (+ heredada legacy_e5_d_...)\n');

%% ========================================================================
%% FIG 7.06  CDF de caudal de los usuarios de BORDE
%%   Regenerada de la figura e5_e de run_e5_adaptive, SIN TITULO.
%% ========================================================================
fprintf('\n--- FIG 7.06  CDF de borde ---\n');
selC = [iR1 iST iAD];
for c = selC
    cdf = E5.KK{c}.edge.cdf_R;
    fprintf('  %-16s CDF con %d puntos | R en [%.4f, %.4f] Mbps\n', ...
        nm5{c}, size(cdf,1), min(cdf(:,1)), max(cdf(:,1)));
end
% El p5 de cada CDF debe coincidir con el R_p5 de borde ya verificado en 7.02
d6 = max(abs(gE(selC) - [5.589765 7.430042 8.267286]));
chk('7.06 R_p5 de borde de las tres curvas', d6, 1e-6);
assert(d6 <= 1e-6, '7.06: las curvas no son las de los casos verificados.');

f = figure('Color','w','Visible','off');
for c = selC
    cdf = E5.KK{c}.edge.cdf_R;
    semilogx(cdf(:,1), cdf(:,2), 'LineWidth',1.7, 'DisplayName',nm5{c}); hold on;
end
yline(0.05,'k:','LineWidth',1.5,'DisplayName','percentil 5');
grid on; box on;
% "caudal", no "throughput": es el termino que usa el resto de la memoria.
xlabel('caudal por usuario R [Mbps]'); ylabel('CDF [-]');
legend('Location','southeast');
bigfont(f); save_fig(f, FIGDIR, 'fig_7_06_cdf_borde'); close(f);
fprintf('  exportada como fig_7_06_cdf_borde\n');

%% ========================================================================
%% FIG 7.A1  Plano (R_agregado, R_p5 borde): el plano DEGENERADO (anexo)
%%   Regenerada de la figura e5_a de run_e5_adaptive, SIN TITULO. Es material de
%%   anexo: en este plano ambas metricas decrecen con alpha, luego alpha->0 domina
%%   y la frontera se reduce a un punto -- la prueba "por encima de la frontera"
%%   es VACUA aqui, que es justamente lo que la figura documenta.
%% ========================================================================
fprintf('\n--- FIG 7.A1  plano degenerado (anexo) ---\n');
nPar = numel(E5.FR.fx);
monoAgg  = all(diff(E5.sw.R_agg) < 0);
monoEdge = all(diff(E5.sw.R_p5_edge) < 0);
fprintf('  puntos no dominados en (R_agg, R_p5_borde): %d de %d (ref 1)\n', nPar, nA);
fprintf('  R_agregado decreciente en alpha: %d | R_p5_borde decreciente: %d (ref 1 y 1)\n', ...
        monoAgg, monoEdge);
assert(nPar == 1 && monoAgg && monoEdge, ...
    '7.A1: el plano ya NO es degenerado; la figura de anexo dejaria de tener sentido.');

f = figure('Color','w','Visible','off');
plot(E5.sw.R_agg/1e3, E5.sw.R_p5_edge, 'o-', 'Color',[.45 .45 .45], 'LineWidth',1.4, ...
     'MarkerFaceColor','w', 'DisplayName','FFR estatica (barrido de \alpha)'); hold on;
for ia = 1:2:nA
    text(E5.sw.R_agg(ia)/1e3, E5.sw.R_p5_edge(ia), sprintf('  %.2f', E5.alphas(ia)), ...
         'FontSize',9, 'Color',[.45 .45 .45]);
end
aggG = cellfun(@(K) K.all.R_agg_Mbps, E5.KK)/1e3;
plot(aggG(iR1), gE(iR1), 'ks', 'MarkerSize',10, 'MarkerFaceColor',[.8 .8 .8], ...
     'DisplayName','reuso-1');
plot(aggG(iST), gE(iST), 'bd', 'MarkerSize',10, 'MarkerFaceColor','b', ...
     'DisplayName',sprintf('FFR estatica (\\alpha*=%.2f)', E5.alpha_star));
plot(aggG(iAD), gE(iAD), 'rp', 'MarkerSize',16, 'MarkerFaceColor','r', ...
     'DisplayName','FFR adaptativa (H2)');
plot(aggG(iIN), gE(iIN), 'v', 'Color',[.85 .5 .1], 'MarkerSize',9, ...
     'MarkerFaceColor',[.85 .5 .1], 'DisplayName','adaptativa INVERTIDA');
plot(E5.FR.fx, E5.FR.fy, 'kp', 'MarkerSize',14, 'MarkerFaceColor','k', ...
     'DisplayName','unico punto no dominado');
grid on; box on;
xlabel('R_{agregado} [Gbps]');
ylabel('R_{p5} de los usuarios de BORDE [Mbps]');
legend('Location','best');
bigfont(f); save_fig(f, FIGDIR, 'fig_7_A1_plano_degenerado'); close(f);
fprintf('  exportada como fig_7_A1_plano_degenerado\n');

%% ---------------------------- cierre -------------------------------------
png = dir(fullfile(FIGDIR,'*.png'));
fprintf('\n=======================================================================\n');
fprintf('  %d PNG en %s\n', numel(png), FIGDIR);
for i = 1:numel(png)
    fprintf('    %-42s %7.1f KB\n', png(i).name, png(i).bytes/1024);
end
fprintf('=======================================================================\n');
fprintf('[cap 7] terminado. Log en %s\n', LOGF);
diary off;

%% ---------------------------- funciones locales --------------------------
function draw_quartiles(Q, qlab, alpha_star, cols)
%DRAW_QUARTILES  Barras de alpha medio por cuartil de elevacion. FUENTE UNICA.
%   Se dibuja en DOS sitios (el panel (a) de fig_7_04 y la figura heredada
%   legacy_fig_7_04_cuartiles) y las barras deben ser IDENTICAS en los dos: tener
%   la regla escrita dos veces es justo el patron que en este proyecto ya causo
%   los bugs de merge_cfg / associate_serving. La leyenda y el ylim los pone cada
%   llamador, que es lo unico que los diferencia.
b = bar(1:4, [Q.nom(:) Q.inv(:) Q.ora(:)], 0.8); hold on;
for k = 1:3, b(k).FaceColor = cols(k,:); end
yline(alpha_star, 'k:', 'LineWidth',1.8);
set(gca,'XTick',1:4,'XTickLabel', qlab);
grid on; box on;
xlabel('cuartil de elevacion del servidor [deg] (instantes ordenados)');
ylabel('\alpha medio, fraccion de banda del INTERIOR [-]');
end
