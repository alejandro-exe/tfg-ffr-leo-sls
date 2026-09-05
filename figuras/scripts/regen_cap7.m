%% REGEN_CAP7   Figuras 7.1 a 7.3 (E5, FFR adaptativa), SIN TITULO.
%
%  Todas CATEGORIA A: lo graficado se LEE de e5_adaptive_results.mat. No se llama a
%  compute_sinr_ffr, compute_interference, compute_kpis, ffr_allocate, ffr_policy,
%  run_one_density ni run_sweep_points.
%
%  Uso: run('.../figuras/scripts/regen_cap7.m')

HERE     = fileparts(mfilename('fullpath'));
PROJROOT = fileparts(fileparts(HERE));
addpath(PROJROOT);
STAGE = fullfile(PROJROOT,'REDACCION','_staging');
FINAL = fullfile(PROJROOT,'REDACCION','figuras');
if ~exist(STAGE,'dir'), mkdir(STAGE); end
promote = @(n) copyfile(fullfile(STAGE,[n '.png']), fullfile(FINAL,[n '.png']));
bigfont = @(fh) set(findall(fh,'-property','FontSize'),'FontSize',11);

E5 = load(fullfile(PROJROOT,'e5_adaptive_results.mat'), ...
          'sw','GR','alphas','KK','cases','alpha_star','DELTA', ...
          'elev_t','alp_t','ctxLite','bin','binName');
iR1 = 1; iST = 2; iAD = 3; iIN = 4;         % run_e5_adaptive.m:180
nm5 = cellfun(@(c) c.name, E5.cases, 'UniformOutput', false);
nA  = numel(E5.alphas);

fprintf('\n################ CAPITULO 7 ################\n');
fprintf('  esquemas: %s | DELTA=%d | alpha* = %.2f\n', strjoin(nm5,' | '), E5.DELTA, E5.alpha_star);

%% ------------- Verificacion numerica contra los valores publicados -------------
% Tabla de E5 (T=1584): R_p5 global / R_p5 BORDE
pub_glob = [5.9 7.4 7.2 7.2];
pub_edge = [5.6 7.4 8.3 8.1];
got_glob = cellfun(@(K) K.all.R_p5_Mbps,  E5.KK);
got_edge = cellfun(@(K) K.edge.R_p5_Mbps, E5.KK);
fprintf('\n--- 7.1 ---\n');
fprintf('  R_p5 global = %s  (publicado %s)\n', mat2str(round(got_glob,2)), mat2str(pub_glob));
fprintf('  R_p5 borde  = %s  (publicado %s)\n', mat2str(round(got_edge,2)), mat2str(pub_edge));
d71a = max(abs(got_glob - pub_glob));
d71b = max(abs(got_edge - pub_edge));
fprintf('  max|dif| global = %.4f | borde = %.4f  (tol 0.05)\n', d71a, d71b);
assert(d71a <= 0.05 && d71b <= 0.05, '7.1 no reproduce la tabla de E5.');
% ganancia de la adaptativa sobre la frontera: publicado -0.063 Mbps
fprintf('  ganancia sobre la frontera = %+.4f Mbps (publicado -0.063)\n', E5.GR.gain_y);
assert(abs(E5.GR.gain_y - (-0.063)) <= 0.001, ...
    '7.1: la ganancia sobre la frontera (%.4f) no coincide con la publicada.', E5.GR.gain_y);
% frontera de Pareto: publicado "11 puntos no dominados de 19"
fprintf('  frontera: %d puntos no dominados de %d (publicado 11 de 19)\n', numel(E5.GR.fx), nA);
assert(numel(E5.GR.fx) == 11 && nA == 19, '7.1: la frontera no tiene 11 de 19 puntos.');
fprintf('  VERIFICACION OK\n');

%% ---------------------------- 7.1 -------------------------------------
f = figure('Color','w','Visible','off');
plot(E5.sw.R_p5, E5.sw.R_p5_edge, 'o-','Color',[.45 .45 .45],'LineWidth',1.4, ...
     'MarkerFaceColor','w','DisplayName','FFR estatica (barrido de \alpha)'); hold on;
plot(E5.GR.fx, E5.GR.fy, 'k-','LineWidth',2.2,'DisplayName','frontera de Pareto estatica');
for ia = 1:2:nA
    text(E5.sw.R_p5(ia), E5.sw.R_p5_edge(ia), sprintf('  %.2f', E5.alphas(ia)), ...
         'FontSize',11,'Color',[.4 .4 .4]);
end
plot(E5.KK{iR1}.all.R_p5_Mbps, E5.KK{iR1}.edge.R_p5_Mbps,'ks','MarkerSize',10, ...
     'MarkerFaceColor',[.8 .8 .8],'DisplayName','reuso-1');
plot(E5.KK{iST}.all.R_p5_Mbps, E5.KK{iST}.edge.R_p5_Mbps,'bd','MarkerSize',10, ...
     'MarkerFaceColor','b','DisplayName',sprintf('FFR estatica (\\alpha*=%.2f)',E5.alpha_star));
plot(E5.KK{iAD}.all.R_p5_Mbps, E5.KK{iAD}.edge.R_p5_Mbps,'rp','MarkerSize',16, ...
     'MarkerFaceColor','r','DisplayName','FFR ADAPTATIVA (H2)');
plot(E5.KK{iIN}.all.R_p5_Mbps, E5.KK{iIN}.edge.R_p5_Mbps,'v','Color',[.85 .5 .1], ...
     'MarkerSize',9,'MarkerFaceColor',[.85 .5 .1],'DisplayName','adaptativa INVERTIDA');
if isfinite(E5.GR.y_frontier)
    plot([E5.GR.x_q E5.GR.x_q],[E5.GR.y_frontier E5.KK{iAD}.edge.R_p5_Mbps],'r--', ...
         'LineWidth',1.2,'DisplayName',sprintf('ganancia = %+.3f Mbps', E5.GR.gain_y));
end
grid on; box on;
xlabel('R_{p5} global [Mbps]'); ylabel('R_{p5} de los usuarios de BORDE [Mbps]');
legend('Location','best');
bigfont(f); save_fig(f, STAGE, 'fig_7_01_pareto_p5global_vs_borde'); close(f);
promote('fig_7_01_pareto_p5global_vs_borde');

%% ---------------------------- 7.2 -------------------------------------
tvec = E5.ctxLite.tvec;
th   = (tvec - tvec(1))/60;
rA   = [min(E5.alp_t)  max(E5.alp_t)];
rE   = [min(E5.elev_t) max(E5.elev_t)];
cc   = corr(E5.alp_t(:), E5.elev_t(:));
fprintf('\n--- 7.2 ---\n');
fprintf('  alpha(t) recorre [%.3f, %.3f]  (publicado [0.468, 0.672])\n', rA);
fprintf('  elev(t)  recorre [%.1f, %.1f] deg (publicado [52.2, 85.5])\n', rE);
fprintf('  corr(alpha, elev) = %+.4f  (publicado +1.000)\n', cc);
assert(max(abs(rA - [0.468 0.672])) <= 0.0006, '7.2: el recorrido de alpha no coincide.');
assert(max(abs(rE - [52.2 85.5]))   <= 0.06,   '7.2: el recorrido de elevacion no coincide.');
assert(abs(cc - 1) <= 5e-4, '7.2: la correlacion alpha-elevacion no es +1.000.');
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off');
yyaxis left
plot(th, E5.elev_t, 'LineWidth',1.6); ylabel('elevacion media del servidor [deg]');
yyaxis right
plot(th, E5.alp_t, 'LineWidth',1.8); hold on;
plot(th, E5.alpha_star*ones(size(th)), '--', 'LineWidth',1.4);
ylabel('\alpha, fraccion de banda del INTERIOR [-]'); ylim([0 1]);
grid on; box on; xlabel('tiempo [min]');
legend({'elevacion','\alpha adaptativa', sprintf('\\alpha estatica = %.2f',E5.alpha_star)}, ...
       'Location','best');
bigfont(f); save_fig(f, STAGE, 'fig_7_02_alpha_vs_elevacion'); close(f);
promote('fig_7_02_alpha_vs_elevacion');

%% ---------------------------- 7.3 -------------------------------------
B = E5.bin.R_p5_edge;                       % [casos x bins]
fprintf('\n--- 7.3 ---\n');
fprintf('  bins: %s\n', strjoin(E5.binName,' | '));
for c = [iR1 iST iAD]
    fprintf('  %-16s %s\n', nm5{c}, mat2str(round(B(c,:),3)));
end
% El tramo 25-35 deg debe estar VACIO a T=1584 (el servidor nunca baja de 49.2 deg)
fprintf('  tramo "%s" vacio: %d (esperado 1)\n', E5.binName{1}, all(isnan(B(:,1))));
assert(all(isnan(B(:,1))), '7.3: el tramo 25-35 deg NO esta vacio a densidad real.');
% adaptativa vs estatica: publicado +15.6 %% en 35-60 y -5.6 %% en 60-90
g2 = 100*(B(iAD,2)/B(iST,2) - 1);
g3 = 100*(B(iAD,3)/B(iST,3) - 1);
fprintf('  adaptativa vs estatica: %s %+.1f %%  |  %s %+.1f %%  (publicado +15.6 / -5.6)\n', ...
        E5.binName{2}, g2, E5.binName{3}, g3);
assert(abs(g2 - 15.6) <= 0.06 && abs(g3 - (-5.6)) <= 0.06, ...
    '7.3: las ganancias por tramo no coinciden con las publicadas.');
fprintf('  VERIFICACION OK\n');

f = figure('Color','w','Visible','off');
bar(B([iR1 iST iAD],:).');
set(gca,'XTickLabel',E5.binName);
grid on; box on;
ylabel('R_{p5} de BORDE [Mbps]'); xlabel('tramo de elevacion del servidor [deg]');
legend(nm5([iR1 iST iAD]),'Location','northwest');
bigfont(f); save_fig(f, STAGE, 'fig_7_03_rp5_por_tramo_elevacion'); close(f);
promote('fig_7_03_rp5_por_tramo_elevacion');

fprintf('\n[cap 7] terminado.\n');
