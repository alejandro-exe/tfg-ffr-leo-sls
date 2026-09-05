%% RUN_E3_BEAMDENSITY  E3: saturacion ESPECTRAL por densificacion de haces.
%
%  MOTIVACION FISICA (el porque de este experimento)
%  --------------------------------------------------
%  La FASE A demostro que la densidad ORBITAL no satura el sistema: los interferentes
%  inter-satelite crecen linealmente con T pero su potencia se queda en I_inter/N <=
%  -21.7 dB incluso con 4000 satelites (el terminal VSAT, ITU-R S.1428, los rechaza
%  fuera de eje), de modo que mas satelites MEJORAN el KPI en vez de degradarlo.
%  La saturacion relevante no es orbital sino ESPECTRAL: el operador estrecha los
%  haces para vender mas capacidad por km2, lo que multiplica las celdas co-canal
%  sobre la MISMA zona y dispara la interferencia INTRA-satelite, que es la
%  dominante segun todos los hallazgos previos (E0, run_ffr_demo, FASE A).
%
%  EJE DE BARRIDO: cfg.radio.beamwidth3dB_deg = [3.0 2.2 1.5 1.0 0.7 0.5] deg.
%  La separacion inter-haz la deriva build_beam_layout con la regla NORMATIVA
%  3GPP TR 38.821 Sec. 6.1.1:  s = rad2deg(sqrt(3)*sin(deg2rad(HPBW)/2)),
%  asi que no hay ningun parametro libre extra: al estrechar el haz, las celdas
%  encogen solas.
%
%  ACOPLAMIENTO FISICO (decision de modelado, cfg.e3.powerMode)
%  -----------------------------------------------------------
%  (1) La ganancia de pico NO puede ser constante al estrechar el haz. Para una
%      apertura circular la directividad va como 1/theta_3dB^2:
%          Gmax(theta) = Gmax_ref + 20*log10(theta_ref/theta)
%      con (theta_ref, Gmax_ref) = (1.5 deg, cfg.radio.Gmax_dBi) el punto actual.
%      OJO (auditoria, hallazgo G2): esa subida de directividad entra en el modelo
%      **a traves de cfg.radio.EIRPdensity_dBWMHz, NO de cfg.radio.Gmax_dBi**. Por
%      definicion la densidad EIRP YA INCLUYE la ganancia de pico del haz, luego
%      sumar Gmax aparte seria doble contabilizacion; de hecho Gmax_dBi se cancela
%      exactamente en compute_link_budget y NO aparece en ffr_context. Gmax(theta)
%      se calcula aqui SOLO para imprimirlo en la tabla del plan (es informativo).
%      El aviso correspondiente esta en la cabecera de compute_link_budget.
%  (2) POTENCIA TOTAL DEL SATELITE CONSTANTE ('total_const', el nominal). Con la
%      misma area de cobertura teselada por celdas mas pequenas hacen falta
%      nBeams_cob ~ 1/theta^2 haces, luego la potencia por haz baja en el MISMO
%      factor en que sube la ganancia:
%          EIRPdens_eff = EIRPdens_ref + (Gmax - Gmax_ref) - 20*log10(theta_ref/theta)
%                       = EIRPdens_ref + 20log10(th_ref/th) - 20log10(th_ref/th)
%                       = EIRPdens_ref                       <-- SE CANCELAN
%      Es decir: en primer orden la subida de directividad y la bajada de potencia
%      por haz se anulan y la PSD de portadora NO cambia con theta. Consecuencia
%      metodologica IMPORTANTE: bajo esta hipotesis el C/N es invariante y TODO lo
%      que se mida al densificar viene de la INTERFERENCIA, que es exactamente lo
%      que se quiere aislar.
%      'per_beam_const' (sensibilidad): la potencia por haz se mantiene y solo sube
%      la ganancia -> EIRPdens_eff = EIRPdens_ref + 20*log10(theta_ref/theta).
%      Es la hipotesis OPTIMISTA (el satelite radia mas potencia total cuanto mas
%      densifica); se incluye para acotar el efecto.                % <-- DECISION
%
%  AVISO METODOLOGICO QUE HAY QUE LEER ANTES DE INTERPRETAR (VARIANTES A y B)
%  --------------------------------------------------------------------------
%  La retícula hexagonal con separacion normativa es AUTO-SIMILAR EN ANGULO: al
%  estrechar el haz, celdas y separacion encogen JUNTAS, de modo que un vecino de
%  primer anillo esta siempre a s/theta ~ 0.866 anchos de haz y aporta siempre
%  -10.67 dB.
%    MATIZ: ese -10.67 es el del vecino ADYACENTE
%    (theta_off = s), primer co-canal SOLO en reuso pleno -- que es el caso de este
%    parrafo ("siempre 90" co-canal = todos los demas haces del cluster). Con Delta=3
%    el primer co-canal esta a sqrt(3)*s (1.500 anchos, -17.954 dB) y con Delta=4 a
%    2*s (1.732 anchos, -18.456 dB). La PROPIEDAD de invariancia se cumple igual de
%    bien en las tres distancias (recorridos 0.010 / 0.014 / 0.028 dB sobre este
%    mismo eje), luego se acota EL NUMERO, no la propiedad: la conclusion de E3 vale
%    para los tres esquemas. El sqrt(3) de u = 1.61634*sqrt(3) viene de la separacion
%    normativa s = sqrt(3)*sin(HPBW/2), NO de la distancia de reuso sqrt(3)*s de
%    Delta=3 (u(s) = 2.79958 frente a u(sqrt(3)*s) = 4.84902).
%  Con nRings FIJO el numero de co-canal tampoco cambia (siempre 90),
%  luego la CIR intra seria practicamente INVARIANTE a theta. Lo que si cambia con
%  nRings fijo es el TAMANO FISICO del cluster (5*s: de ~118 km a theta=3 deg a
%  ~20 km a theta=0.5 deg) frente a la rejilla de usuarios, que es FIJA (40 km):
%  con haces estrechos la rejilla SOBRESALE del cluster y los usuarios de fuera se
%  asocian a un haz muy desapuntado y ven MENOS interferentes de los que les
%  tocarian. Esto ademas es incoherente con la hipotesis (2), que supone el area
%  teselada por ~1/theta^2 haces mientras el modelo de interferencia solo incluye 91.
%
%  Por eso se ejecutan DOS variantes con la MISMA fisica:
%    A) 'nRings=5 FIJO' (la especificada). Se calcula la GUARDA en cada punto y se
%       marca SESGADO lo que no la tenga. NO se silencia.
%    B) 'guarda preservada': nRings(theta) escalado para que el cluster cubra
%       siempre la rejilla + 1 anillo de guarda. Es la variante que DE VERDAD
%       densifica el entorno co-canal y permite separar densificacion de
%       truncamiento del modelo.
%
%  ESQUEMAS (misma poblacion, mismos instantes, misma clasificacion centro/borde):
%    reuse1 | reuseD(3) | reuseD(4) | ffr(3,alpha*) | ffr(4,alpha*) | ffr-adaptativa
%  alpha* = 0.5, el de mejor compromiso del barrido previo (run_ffr_demo: 0.5 en
%  ventana completa; E5: ffr-estatica alpha* = 0.50).
%
%  Uso: ejecutar en la raiz del proyecto. Resultados en e3_results.mat, figuras en
%  figs_e3/. No modifica la fisica del motor (unico cambio de infraestructura:
%  pt.cfg_over en run_one_density). E0 y run_ffr_demo quedan bit-identicos.

clear; clc; close all;

fprintf('=======================================================================\n');
fprintf('  E3: SATURACION ESPECTRAL POR DENSIFICACION DE HACES\n');
fprintf('=======================================================================\n\n');

%% ------------------------------------------------------------------------
%  1. CONFIGURACION BASE (perfil convergido de la FASE A)
%  ------------------------------------------------------------------------
cfg = config_default();

cfg.constellations(1).T = 1584;      % densidad orbital FIJA (no es el eje de E3)
cfg.constellations(1).P = 72;

cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 4;            % FASE A: convergido
cfg.time.dt          = 120;          % FASE A: Nt = 61
cfg.time.duration    = 7200;

cfg.ffr.tau_mode = 'quantile';
cfg.ffr.tau_q    = 50;

cfg.beams.nRings = 5;
cfg.beams.nBeams = [];

cfg.compute.timeBlock = 20;          % el estimador puede pedir bajarlo (ver §3)

% --- Parametros propios de E3 ---
% cfg.e3.powerMode YA VIENE DE config_default ('total_const' por defecto): es una
% hipotesis fisica del sistema, no del experimento, y por eso se movio al perfil.
% Aqui solo se anclan los valores de REFERENCIA del barrido (el punto theta_ref
% respecto del cual se mide dG), que si son del diseno del experimento.
cfg.e3.theta_ref   = cfg.radio.beamwidth3dB_deg;   % 2.08 deg (perfil Ku real)
cfg.e3.Gmax_ref    = cfg.radio.Gmax_dBi;           % 38 dBi (INFORMATIVO, solo se imprime)
cfg.e3.EIRPref     = cfg.radio.EIRPdensity_dBWMHz; % 8.38 dBW/MHz (perfil Ku real)

thetaList = [3.0 2.2 1.5 1.0 0.7 0.5];             % deg, de celdas grandes a densas
alphaStar = 0.5;                                   % del barrido previo (E5)
nWorkers  = [];                                    % [] = automatico

fprintf('Perfil: T=%d P=%d | rejilla %.0f km paso %.0f km | dt=%.0f s | tau cuantil %.0f%%\n', ...
    cfg.constellations(1).T, cfg.constellations(1).P, cfg.ground.radius_km, ...
    cfg.ground.step_km, cfg.time.dt, cfg.ffr.tau_q);
fprintf('Eje: beamwidth3dB = [%s] deg | powerMode = %s | alpha* = %.2f\n\n', ...
    num2str(thetaList), cfg.e3.powerMode, alphaStar);

%% ------------------------------------------------------------------------
%  2. ESQUEMAS
%  ------------------------------------------------------------------------
cases = { ...
    struct('scheme','reuse1','Delta',1,'alpha',NaN,      'adaptive',false,'name','reuse1'), ...
    struct('scheme','reuseD','Delta',3,'alpha',NaN,      'adaptive',false,'name','reuseD(3)'), ...
    struct('scheme','reuseD','Delta',4,'alpha',NaN,      'adaptive',false,'name','reuseD(4)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',alphaStar,'adaptive',false,'name','ffr(3,a*)'), ...
    struct('scheme','ffr',   'Delta',4,'alpha',alphaStar,'adaptive',false,'name','ffr(4,a*)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',alphaStar,'adaptive',true, 'name','ffr-adapt(3)') };
nC = numel(cases);
schemeNames = cellfun(@(c) c.name, cases, 'UniformOutput', false);

%% ------------------------------------------------------------------------
%  3. GEOMETRIA DE CELDAS, ACOPLAMIENTO Y GUARDA POR PUNTO
%  ------------------------------------------------------------------------
% Todo lo que depende SOLO de theta se resuelve aqui, antes de simular, para poder
% imprimir el plan (y las guardas) y decidir si algun punto va a salir sesgado.
nT   = numel(thetaList);
GEO  = struct('theta',[],'s_deg',[],'s_km',[],'cellDiam_km',[],'dens_per1000',[], ...
              'Gmax',[],'EIRP',[],'nRingsA',[],'nRingsB',[],'guardA',[],'guardB',[], ...
              'biasedA',[],'biasedB',[]);
GEO  = repmat(GEO, 1, nT);

h_ref = cfg.constellations(1).h;
Rg    = cfg.ground.radius_km;

fprintf('--------- GEOMETRIA DE CELDAS Y ACOPLAMIENTO ---------\n');
fprintf('%7s %9s %9s %11s %14s %9s %11s %8s %8s %8s\n', ...
    'HPBW','s[deg]','s[km]','diam[km]','haces/1000km2','Gmax','EIRPdens','nR(A)','grd(A)','nR(B)');
for i = 1:nT
    th   = thetaList(i);
    sdeg = rad2deg(sqrt(3) * sin(deg2rad(th)/2));       % 38.821 Sec. 6.1.1
    skm  = h_ref * tand(sdeg);                          % proyeccion nadir (build_beam_layout)

    % Diametro de celda: la celda hexagonal de lado a tiene separacion entre centros
    % s = a*sqrt(3); se reporta el diametro del circulo de igual area.
    Acell   = (sqrt(3)/2) * skm^2;                      % area del hexagono de separacion s
    dcell   = 2*sqrt(Acell/pi);                         % diametro equivalente en area
    dens    = 1000 / Acell;                             % haces por 1000 km2

    % Acoplamiento ganancia / potencia
    dG = 20*log10(cfg.e3.theta_ref / th);               % subida de directividad
    switch lower(cfg.e3.powerMode)
        case 'total_const'
            eirp = cfg.e3.EIRPref + dG - dG;            % SE CANCELAN (explicito a proposito)
        case 'per_beam_const'
            eirp = cfg.e3.EIRPref + dG;
        otherwise
            error('run_e3:powerMode','cfg.e3.powerMode no reconocido: %s', cfg.e3.powerMode);
    end
    % INFORMATIVO SOLAMENTE (auditoria, G2): se imprime en la tabla del plan para
    % documentar la directividad implicada, pero NO se pasa a cfg -- Gmax_dBi no
    % entra en la fisica (esta contenido en EIRPdensity_dBWMHz). El acoplamiento
    % real de este barrido es el `eirp` de arriba.
    gmax = cfg.e3.Gmax_ref + dG;

    % GUARDA. El anillo externo r tiene sus haces mas cercanos a r*s*sqrt(3)/2
    % (puntos medios de lado); se usa esa distancia por ser el caso CONSERVADOR.
    guardOf = @(nR) (nR*skm*sqrt(3)/2 - Rg) / skm;
    gA = guardOf(5);
    % Variante B: menor nRings con guarda >= 1 anillo
    nRB = 2;
    while guardOf(nRB) < 1 && nRB < 40, nRB = nRB + 1; end

    GEO(i).theta = th;   GEO(i).s_deg = sdeg;  GEO(i).s_km = skm;
    GEO(i).cellDiam_km = dcell;  GEO(i).dens_per1000 = dens;
    GEO(i).Gmax = gmax;  GEO(i).EIRP = eirp;
    GEO(i).nRingsA = 5;  GEO(i).guardA = gA;  GEO(i).biasedA = gA <= 0;
    GEO(i).nRingsB = nRB; GEO(i).guardB = guardOf(nRB); GEO(i).biasedB = false;

    fprintf('%7.2f %9.4f %9.2f %11.2f %14.3f %9.2f %11.2f %8d %8.2f %8d\n', ...
        th, sdeg, skm, dcell, dens, gmax, eirp, 5, gA, nRB);
end

nBeamsOf = @(n) 3*n.*(n+1) + 1;
biasedA  = [GEO.biasedA];
fprintf('\nVariante A (nRings=5 FIJO, la especificada): %d de %d puntos SESGADOS por guarda\n', ...
    sum(biasedA), nT);
if any(biasedA)
    fprintf('  SESGADOS: HPBW = [%s] deg -> la rejilla de %.0f km SOBRESALE del cluster.\n', ...
        num2str(thetaList(biasedA)), Rg);
    fprintf('  Esos puntos ven MENOS interferentes de los que les corresponden: su SINR es\n');
    fprintf('  OPTIMISTA y su margen de viabilidad NO es fiable. Se marcan en todas las tablas.\n');
end
fprintf('Variante B (guarda preservada): nRings = [%s] -> nBeams = [%s]\n\n', ...
    num2str([GEO.nRingsB]), num2str(nBeamsOf([GEO.nRingsB])));

%% ------------------------------------------------------------------------
%  4. PUNTOS DEL BARRIDO
%  ------------------------------------------------------------------------
% NOTA (auditoria, G2): el override NO incluye Gmax_dBi. Se sobrescribia antes
% (Gmax_ref + 20log10(th_ref/th)) y era un NO-OP: Gmax_dBi se cancela en
% compute_link_budget y no aparece en ffr_context, porque la ganancia de pico ya va
% dentro de EIRPdensity_dBWMHz. Mantenerlo sugeria un acoplamiento que el codigo no
% implementa; el acoplamiento REAL viaja por EIRPdensity (GEO(i).EIRP).
mkPoint = @(i, nR, tag) struct( ...
    'cfg_over', struct('radio', struct('beamwidth3dB_deg', GEO(i).theta, ...
                                       'EIRPdensity_dBWMHz',GEO(i).EIRP), ...
                       'beams', struct('nRings', nR, 'nBeams', [])), ...
    'cases',       {cases}, ...
    'elMin_valid', 45, ...
    'label',       sprintf('%s-th%.1f', tag, GEO(i).theta));

ptsA = cell(1,nT);  ptsB = cell(1,nT);
for i = 1:nT
    ptsA{i} = mkPoint(i, GEO(i).nRingsA, 'A');
    ptsB{i} = mkPoint(i, GEO(i).nRingsB, 'B');
end

%% ------------------------------------------------------------------------
%  5. PLAN DE MEMORIA (antes de arrancar)
%  ------------------------------------------------------------------------
% El estimador lee nBeams de cfg_base; la variante B usa clusters MAYORES, asi que
% se dimensiona con el PEOR caso para que el plan no se quede corto.
nRworst    = max([GEO.nRingsB]);
cfgWorstOf = @(c) with_rings(c, nRworst);

nCores = max(feature('numcores'), 1);
fprintf('--------- PLAN DE MEMORIA (variante A: %d haces) ---------\n', 3*5*6+1);
MEMA = estimate_sweep_memory(cfg, ptsA, nCores, struct('verbose',true));
fprintf('\n--------- PLAN DE MEMORIA (variante B: peor caso %d haces) ---------\n', ...
    3*max([GEO.nRingsB])*(max([GEO.nRingsB])+1)+1);
MEMB = estimate_sweep_memory(cfgWorstOf(cfg), ptsB, nCores, struct('verbose',true));

% Si con el timeBlock actual no caben los workers, se baja (el troceado es EXACTO:
% test_timeblock_invariance da max|dif| = 0, luego esto no cambia ningun resultado).
tbReco = [MEMA.timeBlock_reco, MEMB.timeBlock_reco];
tbReco = tbReco(~isnan(tbReco));
if ~isempty(tbReco)
    tbNew = min(tbReco);
    fprintf('\n[AVISO] Se reduce cfg.compute.timeBlock de %d a %d para que quepan los workers.\n', ...
        cfg.compute.timeBlock, tbNew);
    cfg.compute.timeBlock = tbNew;
    MEMA = estimate_sweep_memory(cfg, ptsA, nCores, struct('verbose',false));
    MEMB = estimate_sweep_memory(cfgWorstOf(cfg), ptsB, nCores, struct('verbose',false));
end

% El numero de workers lo DECIDE el estimador de memoria corregido (74 B/elem +
% troceado), no el numero de nucleos: en esta maquina el limite es la RAM. Con un
% solo worker se pasa 0 = SECUENCIAL EN EL CLIENTE, que ademas evita reconstruir la
% tabla P.618 dentro de un proceso worker (~90 s) para nada.
nFit = min([MEMA.nFit, MEMB.nFit, nCores]);
if isempty(nWorkers)
    if isnan(nFit) || nFit <= 1
        nWorkers = 0;
    else
        nWorkers = nFit;
    end
end
if nWorkers == 0, modoStr = 'SECUENCIAL (la RAM no da para 2 workers)';
else,             modoStr = sprintf('%d workers', nWorkers); end
fprintf('\n[plan] pico/worker A %.2f GB | B %.2f GB | caben %d -> se ejecuta en %s\n', ...
    MEMA.peak_worker_GB, MEMB.peak_worker_GB, nFit, modoStr);

%% ------------------------------------------------------------------------
%  6. BARRIDO
%  ------------------------------------------------------------------------
fprintf('\n\n########## VARIANTE A: nRings = 5 FIJO (especificada) ##########\n');
[resA, INFOA] = run_sweep_points(cfg, ptsA, nWorkers, struct('verbose',true));

fprintf('\n\n########## VARIANTE B: guarda preservada (control fisico) ##########\n');
[resB, INFOB] = run_sweep_points(cfgWorstOf(cfg), ptsB, nWorkers, struct('verbose',true));

%% ------------------------------------------------------------------------
%  7. TABLAS
%  ------------------------------------------------------------------------
dens = [GEO.dens_per1000];

print_variant_tables('A  (nRings = 5 FIJO)', resA, GEO, cases, dens, biasedA);
print_variant_tables('B  (guarda preservada)', resB, GEO, cases, dens, false(1,nT));

%% ------------------------------------------------------------------------
%  8. MARGEN DE VIABILIDAD (resultado central de H1)
%  ------------------------------------------------------------------------
fprintf('\n\n=======================================================================\n');
fprintf('  MARGEN DE VIABILIDAD: densidad MAXIMA con verdict = ''viable''\n');
fprintf('=======================================================================\n');
fprintf('(la densidad crece al estrechar el haz: el barrido va de %.3f a %.3f haces/1000 km2)\n', ...
    dens(1), dens(end));

MARGIN.A = viability_margin(resA, cases, dens, biasedA, 'A (nRings=5 fijo)');
MARGIN.B = viability_margin(resB, cases, dens, false(1,nT), 'B (guarda preservada)');

%% ------------------------------------------------------------------------
%  9. FIGURAS
%  ------------------------------------------------------------------------
if ~exist('figs_e3','dir'), mkdir('figs_e3'); end
make_e3_figures(resA, resB, cases, dens, biasedA, cfg);

%% ------------------------------------------------------------------------
%  10. Guardado
%  ------------------------------------------------------------------------
save('e3_results.mat', 'resA', 'resB', 'INFOA', 'INFOB', 'GEO', 'cases', ...
     'schemeNames', 'thetaList', 'dens', 'biasedA', 'cfg', 'alphaStar', ...
     'MARGIN', 'MEMA', 'MEMB', '-v7.3');
fprintf('\nResultados en e3_results.mat | figuras en figs_e3/\n');

%% ========================================================================
%  FUNCIONES LOCALES
%  ========================================================================
function c = with_rings(c, nR)
%WITH_RINGS  Copia de cfg con el cluster fijado por anillos completos.
c.beams.nRings = nR;
c.beams.nBeams = [];
end

% -------------------------------------------------------------------------
function print_variant_tables(title, res, GEO, cases, dens, biased)
%PRINT_VARIANT_TABLES  Tablas por esquema de una variante del barrido.
nT = numel(res);  nC = numel(cases);
fprintf('\n\n=======================================================================\n');
fprintf('  VARIANTE %s\n', title);
fprintf('=======================================================================\n');

for c = 1:nC
    fprintf('\n--- %s ---\n', cases{c}.name);
    fprintf('%7s %9s %8s %10s %9s %9s %9s %9s %8s %8s %9s %10s %s\n', ...
        'HPBW','dens','nBeams','SINRe_p5','cobert','margFlr','Rp5_bor','Rp5_glob', ...
        'R_agg','Jain','SE_e_p5','penaliz','veredicto');
    for i = 1:nT
        K = res{i}.K{c};  v = K.viab;
        if biased(i), mk = '  SESGADO'; else, mk = ''; end
        fprintf('%7.2f %9.3f %8d %10.3f %9.2f %9.3f %9.2f %9.2f %8.2f %8.3f %9.3f %10.2f %-9s%s\n', ...
            GEO(i).theta, dens(i), res{i}.nBeams, v.SINR_edge_p5, v.coverage, ...
            v.margin_floor_dB, K.edge.R_p5_Mbps, K.all.R_p5_Mbps, ...
            K.all.R_agg_Mbps/1e3, K.all.jain, v.SE_edge_p5, K.penalty_mean, ...
            v.verdict, mk);
    end
end

% Desglose intra/inter de la penalizacion (diagnostico del regimen)
fprintf('\n--- Penalizacion por interferencia (C/N - SINR) [dB], media ---\n');
fprintf('%7s %9s', 'HPBW', 'dens');
for c = 1:nC, fprintf(' %13s', cases{c}.name); end
fprintf('\n');
for i = 1:nT
    fprintf('%7.2f %9.3f', GEO(i).theta, dens(i));
    for c = 1:nC, fprintf(' %13.2f', res{i}.K{c}.penalty_mean); end
    if biased(i), fprintf('   SESGADO'); end
    fprintf('\n');
end
end

% -------------------------------------------------------------------------
function M = viability_margin(res, cases, dens, biased, tag)
%VIABILITY_MARGIN  Densidad maxima con verdict='viable' por esquema.
%   Devuelve tambien la densidad de CRUCE interpolada con margin_th_dB
%   (= SINR_edge_p5 - gamma_th) en log10(densidad), que da una cifra continua en
%   vez de saltar entre puntos del barrido.
nT = numel(res);  nC = numel(cases);
fprintf('\n--- Variante %s ---\n', tag);
fprintf('%-14s %14s %16s %14s %s\n', 'esquema','dens max viable','dens cruce interp','margen 0dB','veredictos (de menos a mas denso)');

M = struct('scheme',{},'dens_max_viable',{},'dens_cross',{},'verdicts',{});
for c = 1:nC
    ver = cell(1,nT);  mth = nan(1,nT);
    for i = 1:nT
        ver{i} = res{i}.K{c}.viab.verdict;
        mth(i) = res{i}.K{c}.viab.margin_th_dB;
    end
    okv = strcmp(ver,'viable') & ~biased;     % los sesgados NO cuentan
    if any(okv), dmax = max(dens(okv)); else, dmax = NaN; end

    % Cruce interpolado de margin_th_dB = 0 sobre los puntos NO sesgados
    dcross = NaN;
    ii = find(~biased);
    for k = 1:numel(ii)-1
        a = ii(k); b = ii(k+1);
        if ~isnan(mth(a)) && ~isnan(mth(b)) && (mth(a) > 0) && (mth(b) <= 0)
            la = log10(dens(a)); lb = log10(dens(b));
            t  = mth(a) / (mth(a) - mth(b));
            dcross = 10^(la + t*(lb - la));
            break;
        end
    end

    vs = strjoin(cellfun(@(s) s(1:3), ver, 'UniformOutput', false), '>');
    fprintf('%-14s %14s %16s %14.3f %s\n', cases{c}.name, ...
        num2str(dmax,'%.3f'), num2str(dcross,'%.3f'), mth(end), vs);

    M(c).scheme = cases{c}.name;
    M(c).dens_max_viable = dmax;
    M(c).dens_cross = dcross;
    M(c).verdicts = ver;
end

% Ganancia de la FFR sobre reuse1 en el margen
i1 = find(strcmp({M.scheme},'reuse1'),1);
for c = 1:nC
    if c == i1, continue; end
    if ~isnan(M(c).dens_cross) && ~isnan(M(i1).dens_cross)
        fprintf('  %s vs reuse1: margen x%.2f en densidad de cruce\n', ...
            M(c).scheme, M(c).dens_cross / M(i1).dens_cross);
    end
end
end

% -------------------------------------------------------------------------
function make_e3_figures(resA, resB, cases, dens, biasedA, cfg)
%MAKE_E3_FIGURES  Las cuatro figuras pedidas para E3.
nT = numel(resA);  nC = numel(cases);
col = lines(nC);
gth = cfg.viab.gamma_th_dB;  gfl = cfg.viab.gamma_floor_dB;
nm  = cellfun(@(c) c.name, cases, 'UniformOutput', false);

getf_ = @(res, c, f) arrayfun(@(i) subsref_kpi(res{i}.K{c}, f), 1:nT);

% (a) SINR_edge_p5 vs densidad, con los umbrales del criterio de viabilidad
figure('Name','E3(a) SINR de borde vs densidad de haces','Position',[80 80 1100 450]);
for v = 1:2
    subplot(1,2,v);
    if v == 1, res = resA; ttl = 'A: nRings = 5 FIJO'; else, res = resB; ttl = 'B: guarda preservada'; end
    hold on; grid on;
    for c = 1:nC
        plot(dens, getf_(res,c,'viab.SINR_edge_p5'), '-o', 'Color', col(c,:), ...
            'LineWidth', 1.5, 'MarkerFaceColor', col(c,:), 'MarkerSize', 4);
    end
    yline(gth, 'k--', sprintf('\\gamma_{th} = %.1f dB', gth), 'LineWidth', 1.2);
    yline(gfl, 'r--', sprintf('\\gamma_{floor} = %.1f dB', gfl), 'LineWidth', 1.2);
    if v == 1
        for i = find(biasedA), xline(dens(i), ':', 'Color', [.6 .6 .6]); end
    end
    set(gca,'XScale','log');
    xlabel('densidad de haces [haces / 1000 km^2]'); ylabel('SINR_{edge,p5} [dB]');
    title(ttl); if v == 1, legend(nm, 'Location','southwest','FontSize',7); end
end
sgtitle('E3(a): SINR de borde vs densidad de haces (lineas grises = puntos SESGADOS)');
save_fig(gcf, 'figs_e3', 'e3_a_sinr_edge_vs_density');

% (b) Mapa de veredictos
figure('Name','E3(b) Veredicto vs densidad','Position',[80 80 1100 450]);
lev = containers.Map({'inviable','marginal','viable','sin-datos'}, {1,2,3,0});
for v = 1:2
    subplot(1,2,v);
    if v == 1, res = resA; ttl = 'A: nRings = 5 FIJO'; else, res = resB; ttl = 'B: guarda preservada'; end
    Mv = zeros(nC, nT);
    for c = 1:nC, for i = 1:nT, Mv(c,i) = lev(res{i}.K{c}.viab.verdict); end, end
    imagesc(Mv, [0 3]); colormap([0.85 0.85 0.85; 0.80 0.20 0.20; 0.95 0.75 0.20; 0.20 0.65 0.30]);
    set(gca,'YTick',1:nC,'YTickLabel',nm,'FontSize',8);
    set(gca,'XTick',1:nT,'XTickLabel',arrayfun(@(d) sprintf('%.2f',d), dens, 'UniformOutput', false));
    xlabel('densidad [haces / 1000 km^2]'); title(ttl);
    if v == 1
        for i = find(biasedA)
            text(i, 0.35, 'SESG', 'HorizontalAlignment','center','FontSize',7,'Color',[.3 .3 .3]);
        end
    end
end
sgtitle('E3(b): VEREDICTO de viabilidad (rojo=inviable, ambar=marginal, verde=viable)');
save_fig(gcf, 'figs_e3', 'e3_b_verdict_map');

% (c) Compromiso: R_p5 de borde y R agregado
figure('Name','E3(c) Compromiso borde/agregado','Position',[80 80 1100 450]);
subplot(1,2,1); hold on; grid on;
for c = 1:nC
    plot(dens, getf_(resB,c,'edge.R_p5_Mbps'), '-o','Color',col(c,:),'LineWidth',1.5,'MarkerSize',4);
end
set(gca,'XScale','log'); xlabel('densidad [haces / 1000 km^2]'); ylabel('R_{p5} de BORDE [Mbps]');
title('Throughput de borde (variante B)'); legend(nm,'Location','best','FontSize',7);
subplot(1,2,2); hold on; grid on;
for c = 1:nC
    plot(dens, getf_(resB,c,'all.R_agg_Mbps')/1e3, '-o','Color',col(c,:),'LineWidth',1.5,'MarkerSize',4);
end
set(gca,'XScale','log'); xlabel('densidad [haces / 1000 km^2]'); ylabel('R agregado [Gbps]');
title('Caudal agregado (variante B)');
sgtitle('E3(c): compromiso borde vs agregado al densificar');
save_fig(gcf, 'figs_e3', 'e3_c_tradeoff');

% (d) Penalizacion por interferencia
figure('Name','E3(d) Penalizacion por interferencia','Position',[80 80 1100 450]);
for v = 1:2
    subplot(1,2,v);
    if v == 1, res = resA; ttl = 'A: nRings = 5 FIJO'; else, res = resB; ttl = 'B: guarda preservada'; end
    hold on; grid on;
    for c = 1:nC
        plot(dens, arrayfun(@(i) res{i}.K{c}.penalty_mean, 1:nT), '-o', ...
            'Color', col(c,:), 'LineWidth', 1.5, 'MarkerSize', 4);
    end
    set(gca,'XScale','log');
    xlabel('densidad [haces / 1000 km^2]'); ylabel('C/N - SINR [dB]');
    title(ttl); if v == 1, legend(nm,'Location','best','FontSize',7); end
end
sgtitle('E3(d): penalizacion por interferencia vs densidad (debe CRECER si densificar satura)');
save_fig(gcf, 'figs_e3', 'e3_d_penalty');

fprintf('\n4 figuras exportadas a figs_e3/\n');
end

% -------------------------------------------------------------------------
function v = subsref_kpi(K, f)
%SUBSREF_KPI  Lee un campo anidado 'a.b' de la estructura de KPIs.
p = strsplit(f, '.');
v = K;
for i = 1:numel(p), v = v.(p{i}); end
end
