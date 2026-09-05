%% AUDITAR_REGENERABILIDAD   TAREA 1 del encargo 2
%  Para cada figura del inventario, contrasta UNA A UNA las variables que su bloque
%  de graficado necesita contra las que su .mat contiene realmente, y la clasifica:
%
%    A - regenerable directa            : el .mat contiene todo.
%    B - regenerable con re-derivacion  : falta algo, pero es puramente GEOMETRICO
%                                         y reconstruible desde cfg / elementos guardados.
%    C - no regenerable                 : falta algo que exigiria RECALCULAR resultados,
%                                         o no hay .mat.
%
%  Solo LEE. No escribe ficheros, no grafica, no simula.
%
%  Uso: run('.../figuras/scripts/auditar_regenerabilidad.m')

HERE     = fileparts(mfilename('fullpath'));
PROJROOT = fileparts(fileparts(HERE));
addpath(PROJROOT);

% {figura, mat, {especificaciones que la figura NECESITA}, {lo que NO esta y como se obtiene}}
INV = {
 '5.1  calib vs 3GPP',      'calibration_e0_results.mat', ...
    {'S.CIR.stats.CIR_p5','S.CIR.stats.CIR_p50','S.CIR.stats.CIR_p95'}
 % 5.2: lo que la figura dibuja es V = validate_geometry_satscenario(cfg), que NO
 % se guarda en el .mat. Se listan sus campos a proposito para que salgan FALTA.
 '5.2  validacion geom',    'calibration_e0_results.mat', ...
    {'S.cfg.time.dt','S.V.el_own','S.V.el_tb','S.V.rng_own','S.V.rng_tb','S.V.jsat'}
 '5.3  converg. rejilla',   'convergence_studyA.mat', ...
    {'S.ST.steps_km','S.T1.ok','S.T1.p5'}
 '5.4  suelo vs Bessel',    'convergence_nrings.mat', ...
    {'S.RES{end}.spacing_deg','S.RES{end}.nRings','S.cfg.radio.beamwidth3dB_deg', ...
     'S.cfg.radio.sidelobe_floor_dB'}
 '6.1  CDF thr. borde',     'ffr_results.mat', ...
    {'S.KK{1}.edge.cdf_R','S.KK{5}.edge.cdf_R','S.cases{1}.name'}
 '6.2  CDF SINR',           'ffr_results.mat', ...
    {'S.KK{1}.all.cdf_SINR','S.cfg.kpi.gamma0_dB','S.cases{1}.name'}
 '6.3  interf. vs T',       'convergence_studyA.mat', ...
    {'S.ST.T_list','S.T4.ok','S.T4.nCo','S.T4.p5'}
 '6.4  SINR vs densidad',   'e3_results.mat', ...
    {'S.dens','S.resA{1}.K{1}.viab.SINR_edge_p5','S.resB{1}.K{1}.viab.SINR_edge_p5', ...
     'S.biasedA','S.cases{1}.name','S.cfg.viab.gamma_th_dB','S.cfg.viab.gamma_floor_dB'}
 '6.5  mapa veredictos T',  'e3_results.mat', ...
    {'S.dens','S.resA{1}.K{1}.viab.verdict','S.resB{1}.K{1}.viab.verdict', ...
     'S.biasedA','S.cases{1}.name'}
 '6.6  SINR vs mascara',    'minelev_results.mat', ...
    {'S.minElevList','S.R{1}{1}.K{1}.viab.SINR_edge_p5','S.densList(1).name', ...
     'S.densList(1).T','S.nm{1}','S.cfg.viab.gamma_th_dB'}
 '6.7  mapa veredictos me', 'minelev_results.mat', ...
    {'S.minElevList','S.R{1}{1}.K{1}.viab.verdict','S.densList(1).name','S.nm{1}'}
 '6.8  KPI vs operadores',  'interconstellation_results.mat', ...
    {'S.nOpsList','S.res{1}.K{1}.viab.SINR_edge_p1','S.res{1}.K{1}.viab.SINR_edge_p5', ...
     'S.nm{1}','S.cfg.viab.gamma_th_dB'}
 '7.1  frontera Pareto',    'e5_adaptive_results.mat', ...
    {'S.sw.R_p5','S.sw.R_p5_edge','S.GR.fx','S.GR.fy','S.GR.y_frontier','S.GR.x_q', ...
     'S.GR.gain_y','S.alphas','S.KK{1}.all.R_p5_Mbps','S.KK{1}.edge.R_p5_Mbps', ...
     'S.alpha_star','S.DELTA','S.cases{1}.name'}
 '7.2  alpha vs elevacion', 'e5_adaptive_results.mat', ...
    {'S.elev_t','S.alp_t','S.alpha_star','S.ctxLite.tvec'}
 '7.3  Rp5 por tramo',      'e5_adaptive_results.mat', ...
    {'S.bin.R_p5_edge','S.binName','S.cases{1}.name'}
 '8.1  compresion angular', 'minelev_results.mat', ...
    {'S.elGrid','S.sepRad','S.sepTan','S.sepNom','S.STRAT.bins','S.STRAT.compr','S.STRAT.SINR'}
 '8.2  autosimilitud',      'oneweb_results.mat', ...
    {'S.AS.HPBW','S.AS.dev','S.cfg.radio.beamwidth3dB_deg'}
 '8.3  composicion p5',     'ffr_results.mat', ...
    {'S.p5split(1).fracC_R','S.KK{1}.centerFrac','S.cases{1}.name'}
 '8.4  topologia hex/rect', 'oneweb_results.mat', ...
    {'S.OUT.K{1}.viab.SINR_edge_p5','S.OUTR.K{1}.viab.SINR_edge_p5','S.nm{1}'}
 '8.5a oneweb elevacion',   'oneweb_results.mat', ...
    {'S.ELIST','S.RE{1}.K{1}.viab.SINR_edge_p5','S.nm{1}'}
 '8.5b oneweb densific.',   'oneweb_results.mat', ...
    {'S.densB','S.RD{1}.K{1}.viab.SINR_edge_p5','S.nm{1}'}
};

fprintf('\n======================= TAREA 1: REGENERABILIDAD =======================\n');

matCache = containers.Map();
for k = 1:size(INV,1)
    fig = INV{k,1};  matf = INV{k,2};  need = INV{k,3};

    if ~isKey(matCache, matf)
        p = fullfile(PROJROOT, matf);
        if ~exist(p,'file')
            fprintf('\n--- %s ---\n  .mat: %s  ** NO EXISTE ** -> categoria C\n', fig, matf);
            continue
        end
        matCache(matf) = load(p);
    end
    S = matCache(matf); %#ok<NASGU>  (lo usan las cadenas de `need` via eval)

    fprintf('\n--- %s   [%s] ---\n', fig, matf);
    nOK = 0;  falta = {};
    for j = 1:numel(need)
        try
            v = eval(need{j});
            if isempty(v) && ~ischar(v)
                fprintf('    %-46s VACIO\n', need{j});
                falta{end+1} = need{j}; %#ok<SAGROW>
            else
                fprintf('    %-46s OK   %s %s\n', need{j}, class(v), mat2str(size(v)));
                nOK = nOK + 1;
            end
        catch
            fprintf('    %-46s ** FALTA **\n', need{j});
            falta{end+1} = need{j}; %#ok<SAGROW>
        end
    end
    if isempty(falta)
        fprintf('  => %d/%d presentes. CATEGORIA A (regenerable directa)\n', nOK, numel(need));
    else
        fprintf('  => faltan %d: %s\n', numel(falta), strjoin(falta, ', '));
    end
end

fprintf('\n=======================================================================\n');
fprintf(['NOTA: las constantes de presentacion (col = lines(nC), nC = numel(cases),\n' ...
         'nm = cellfun(@(c)c.name,cases), gth/gfl = cfg.viab.*, indices iR1..iIN = 1..4)\n' ...
         'no se auditan: son deterministas y no son datos.\n']);
