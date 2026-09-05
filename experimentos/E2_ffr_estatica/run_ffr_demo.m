%% RUN_FFR_DEMO  FFR sobre el layout multihaz Earth-fixed: reuso-1 vs reuso-D vs FFR.
%  RUNNER DELGADO del perfil de ESTUDIO (Ku, config_default). Toda la fisica esta
%  en el MOTOR y se comparte con E0 y con los futuros experimentos:
%
%     build_beam_layout -> ffr_context -> ffr_coloring -> ffr_policy
%                       -> ffr_allocate -> compute_sinr_ffr -> compute_kpis
%
%  Compara sobre la MISMA poblacion de usuarios e instantes:
%     reuse1                  (peor caso espectral: todos los haces co-canal)
%     reuseD  Delta = 3 y 4   (proteccion maxima: solo co-canal el mismo color)
%     ffr     Delta = 3 y 4   (alpha = cfg.ffr.alpha: interior en banda plena,
%                              borde en su sub-banda de color)
%  y hace un barrido de alpha 0.1:0.1:0.9 para trazar la CURVA DE COMPROMISO
%  (R_agregado vs R_p5 de borde).
%
%  VALIDACION esperada (H1): al pasar de reuse1 a ffr, el R_p5 de los usuarios de
%  BORDE debe SALTAR (es la metrica reina) a costa de ancho de banda; reuseD da la
%  proteccion maxima y ffr debe quedar en un punto intermedio MEJOR en compromiso.
%
%  Uso: ejecutar en la raiz del proyecto.  Resultados en ffr_results.mat.

clear; clc; close all;

% CRONOMETRO DE EXTREMO A EXTREMO. El warmup de la
% tabla P.618 se mide APARTE y se descuenta: es un coste fijo de la configuracion
% atmosferica (~80 s), esta cacheado en un 'persistent' y no depende del
% experimento, luego incluirlo en el tiempo del runner falsearia la comparacion
% entre experimentos.
tRun  = tic;
FIGDIR = 'figs_e2';                 % E2 = comparativa de esquemas (run_ffr_demo)

%% 1. Perfil de estudio: rejilla LOCAL + constelacion reducida (rapidez)
cfg = config_default();

% CONSTELACION REAL Starlink Shell-1 (Walker 53:1584/72/1), la de config_default.
% ANTES este runner usaba T=66/P=6 con el argumento de que "la interferencia
% relevante es la INTRA-satelite, asi que no hace falta saturar". El argumento es
% correcto en cuanto al MECANISMO (y sigue siendo el hallazgo de FASE A: la
% interferencia inter-satelite es despreciable, I_inter/N <= -26 dB incluso con
% 4000 satelites), pero NO justifica publicar las cifras de cabecera de H1 sobre
% una constelacion de juguete: a T=66 la cobertura es del 44% y el servidor pasa
% mucho tiempo a elevacion baja, lo que obligaba a filtrar por "geometria valida"
% (elev >= 45 deg) para dar una tabla defendible. A densidad REAL la cobertura es
% ~1 y el servidor nunca baja de ~52 deg (hallazgo de E3b), luego el filtro deja
% de hacer falta y la tabla de la ventana completa YA es la tabla defendible.
cfg.constellations(1).T = 1584;
cfg.constellations(1).P = 72;

% Rejilla local de usuarios (centro y borde de celda)
cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 3;

% dt grande (el coste es M*N*Nt, ver aviso de build_user_grid) y ventana de 2 h.
% La ventana LARGA es necesaria: con T=66 y 30 min el servidor no pasa nunca de
% 31 deg de elevacion y el cluster Earth-fixed sale muy comprimido angularmente
% (interferencia patologica). Con 2 h se recorren elevaciones de 25 a ~77 deg.
cfg.time.dt       = 60;        % s
cfg.time.duration = 7200;      % s (2 h)

% Umbral centro/borde por CUANTIL (ver ffr_policy): el nivel ABSOLUTO de SINR en
% reuso-1 depende fuertisimamente de la elevacion, asi que un tau fijo (p.ej. los
% 5 dB de config_default) deja el 100% de usuarios en un lado y la FFR DEGENERA
% (todos borde -> peor que reuso-Delta). Fijamos el REPARTO (50/50) y no el nivel.
% El valor absoluto de tau queda ligado al umbral de viabilidad de H1 (tarea abierta).
cfg.ffr.tau_mode = 'quantile';
cfg.ffr.tau_q    = 50;         % % de usuarios de CENTRO

% Cluster de haces: 5 anillos completos = 91 haces. Se elige asi para que la
% rejilla de usuarios (40 km ~ 3.2 celdas) quede con >= 1.5 anillos de GUARDA de
% celdas co-canal alrededor; si no, los usuarios del borde de la rejilla verian
% menos interferentes de los que les corresponden y la FFR pareceria mejor.
cfg.beams.nRings = 5;

% TROCEADO TEMPORAL: OBLIGATORIO a densidad real. El pico del pipeline es
% [M x N x Ntb] dentro de compute_interference (~74 B/elem, mem_peak_model_GB):
% con M=553, Nt=121, N=1584 y 91 haces son 7.53 GB de una pieza, MAS que la RAM
% total de la maquina de referencia (7.45 GB). El troceado es EXACTO (el pipeline
% no acopla instantes; verificado en test_timeblock_invariance y, aqui, contra el
% .mat publicado a T=66 con max|dif| = 0), luego no cambia ningun resultado: solo
% acota el pico. El valor se fija abajo con estimate_sweep_memory.
cfg.compute.timeBlock = 10;                                        % <-- DECISION

fprintf('===== FFR DEMO (perfil de ESTUDIO Ku) =====\n');
fprintf('Ku DL: f=%.0f GHz | B=%.0f MHz | EIRPdens=%.1f dBW/MHz | HPBW=%.2f deg | G/T=%.1f dB/K\n', ...
    cfg.radio.freq_GHz, cfg.radio.B_MHz, cfg.radio.EIRPdensity_dBWMHz, ...
    cfg.radio.beamwidth3dB_deg, cfg.radio.GT_dBK);

% Warmup EXPLICITO de la tabla P.618 (coste fijo, fuera del cronometro del
% experimento). No cambia ningun numero: test_atm_cache verifica que la version
% cacheada es identica bit a bit a la no cacheada.
tW = tic;  atm_loss_dB(cfg, 45);  tWarm = toc(tW);
fprintf('[warmup] tabla P.618 construida en %.1f s (coste fijo, se descuenta)\n', tWarm);
tExp = tic;                        % <- cronometro del EXPERIMENTO, ya con cache

%% 2. DIMENSIONADO DE MEMORIA (antes de simular nada)
%  estimate_sweep_memory es el mismo estimador que dimensiona los barridos (fuente
%  unica del modelo en mem_peak_model_GB). Aqui el "barrido" es de un solo punto y
%  un solo proceso, asi que lo que interesa es .peak_GB y el .timeBlock_reco.
MEM = estimate_sweep_memory(cfg, {struct('T', cfg.constellations(1).T)}, 1, ...
                            struct('verbose', true));
if ~isempty(MEM.timeBlock_reco) && ~isnan(MEM.timeBlock_reco)
    fprintf(['[memoria] timeBlock recomendado por el estimador: %d ' ...
             '(configurado: %d)\n'], MEM.timeBlock_reco, cfg.compute.timeBlock);
end

%% 3. Pipeline por BLOQUES TEMPORALES -> contexto compartido del motor
%  build_ctx_blocked recorre la ventana en bloques ejecutando el MISMO pipeline
%  (propagate -> geometria -> servidor -> enlace -> interferencia -> ffr_context) y
%  devuelve el ctx COMPLETO, identico al de la ventana de una pieza. Todo lo que
%  viene despues (esquemas, barrido de alpha, KPIs, figuras) no se entera.
[ctx, sats, users, BL, BINFO] = build_ctx_blocked(cfg, true);

fprintf('Rejilla local: %d usuarios (radio %.0f km, paso %.0f km) | Nt=%d | N=%d satelites\n', ...
    users.M, cfg.ground.radius_km, cfg.ground.step_km, ctx.Nt, sats.N);
fprintf('\nLayout: %d haces | s = %.4f deg = %.2f km (3GPP TR 38.821 Sec.6.1.1)\n', ...
    BL.nBeams, BL.spacing_deg, BL.spacing_km);

tvec = ctx.tvec;

%% 4. Esquemas a comparar
cases = { ...
    struct('scheme','reuse1','Delta',1,'alpha',NaN, 'name','reuse1'), ...
    struct('scheme','reuseD','Delta',3,'alpha',NaN, 'name','reuseD (D=3)'), ...
    struct('scheme','reuseD','Delta',4,'alpha',NaN, 'name','reuseD (D=4)'), ...
    struct('scheme','ffr',   'Delta',3,'alpha',0.4, 'name','FFR (D=3, a=0.4)'), ...
    struct('scheme','ffr',   'Delta',4,'alpha',0.4, 'name','FFR (D=4, a=0.4)') };
nC = numel(cases);
KK = cell(1,nC);  AA = cell(1,nC);  FF = cell(1,nC);

fprintf('\n--------- EVALUACION POR ESQUEMA ---------\n');
for c = 1:nC
    cc = cases{c};
    cfg.ffr.scheme = cc.scheme;
    cfg.ffr.Delta  = cc.Delta;
    if ~isnan(cc.alpha), cfg.ffr.alpha = cc.alpha; end

    fprintf('\n>> %s\n', cc.name);
    AA{c} = ffr_allocate(cfg, ctx);
    FF{c} = compute_sinr_ffr(cfg, ctx, AA{c});
    KK{c} = compute_kpis(cfg, FF{c}, AA{c});

    % La SINR de reuso-1 es la REFERENCIA del clasificador 'sinr': se calcula una
    % sola vez y se guarda en ctx para que los demas esquemas la reutilicen (y
    % para que la clasificacion centro/borde sea IDENTICA en todos ellos).
    if isempty(ctx.SINR_ref_dB)
        ctx.SINR_ref_dB = AA{c}.SINR_ref_dB;
    end
end

%% 4b. COMPROBACION: la clasificacion CENTRO/BORDE es COMUN a los tres esquemas
%  Es lo que hace justa la comparacion del criterio de viabilidad: reuse1 y
%  reuseD no tienen particion centro/borde nativa, asi que se les aplica la MISMA
%  clasificacion de referencia que a la FFR (SINR de reuso-1 frente a tau). No es
%  una decision de este runner: ffr_allocate clasifica ANTES de ramificar por
%  esquema, luego isCenter no depende de cfg.ffr.scheme. Aqui se VERIFICA.
for c = 2:nC
    if ~isequal(AA{c}.isCenter, AA{1}.isCenter)
        error('run_ffr_demo:edgePop', ...
            ['La clasificacion centro/borde NO coincide entre "%s" y "%s": los KPIs ' ...
             'de borde y el veredicto de viabilidad compararian poblaciones distintas.'], ...
            cases{1}.name, cases{c}.name);
    end
end
fprintf(['\n[check] clasificacion centro/borde IDENTICA en los %d esquemas ' ...
         '(%d muestras de borde de %d con cobertura): la viabilidad se juzga sobre ' ...
         'la MISMA poblacion.\n'], nC, sum(ctx.cov(:) & ~AA{1}.isCenter(:)), sum(ctx.cov(:)));

%% 5. (c) Tabla comparativa
fprintf('\n================================ TABLA COMPARATIVA ================================\n');
fprintf('%-18s %7s %8s %9s %9s %9s %7s %7s\n', ...
    'esquema','B[MHz]','SINR_p5','R_p5','R_media','R_agreg','P(cob)','Jain');
fprintf('%-18s %7s %8s %9s %9s %9s %7s %7s\n', ...
    '','','[dB]','[Mbps]','[Mbps]','[Gbps]','','');
fprintf('%s\n', repmat('-',1,82));
for c = 1:nC
    K = KK{c};
    fprintf('%-18s %7.1f %8.2f %9.1f %9.1f %9.2f %7.2f %7.3f\n', ...
        cases{c}.name, K.B_user_MHz_mean, K.all.SINR_p5, K.all.R_p5_Mbps, ...
        K.all.R_mean_Mbps, K.all.R_agg_Mbps/1e3, K.all.Pcov, K.all.jain);
end
fprintf('%s\n', repmat('-',1,82));
fprintf('DESGLOSE CENTRO / BORDE (clasificador ''%s'', umbral %.2f):\n', ...
    cfg.ffr.classifier, cfg.ffr.tau_dB);
fprintf('%-18s %9s %9s %9s %9s %9s\n', ...
    'esquema','SINRp5_C','SINRp5_B','R_p5_C','R_p5_B','%centro');
fprintf('%s\n', repmat('-',1,68));
for c = 1:nC
    K = KK{c};
    fprintf('%-18s %9.2f %9.2f %9.1f %9.1f %9.1f\n', ...
        cases{c}.name, K.center.SINR_p5, K.edge.SINR_p5, ...
        K.center.R_p5_Mbps, K.edge.R_p5_Mbps, 100*K.centerFrac);
end
fprintf('%s\n', repmat('-',1,68));

% Salto del R_p5 de BORDE de reuse1 a FFR (resultado clave de H1)
iR1  = find(strcmpi(cellfun(@(x)x.scheme,cases,'UniformOutput',false),'reuse1'), 1);
iF3  = find(cellfun(@(x)strcmpi(x.scheme,'ffr') && x.Delta==3, cases), 1);
iF4  = find(cellfun(@(x)strcmpi(x.scheme,'ffr') && x.Delta==4, cases), 1);
e1   = KK{iR1}.edge.R_p5_Mbps;
fprintf('\n[H1] R_p5 de BORDE: reuse1 = %.1f Mbps -> FFR(D=3) = %.1f Mbps (x%.2f, %+.1f Mbps)\n', ...
    e1, KK{iF3}.edge.R_p5_Mbps, KK{iF3}.edge.R_p5_Mbps/max(e1,eps), KK{iF3}.edge.R_p5_Mbps-e1);
fprintf('[H1] R_p5 de BORDE: reuse1 = %.1f Mbps -> FFR(D=4) = %.1f Mbps (x%.2f, %+.1f Mbps)\n', ...
    e1, KK{iF4}.edge.R_p5_Mbps, KK{iF4}.edge.R_p5_Mbps/max(e1,eps), KK{iF4}.edge.R_p5_Mbps-e1);
fprintf('[H1] SINR_p5 de BORDE: reuse1 = %.2f dB -> FFR(D=3) = %.2f dB (%+.2f dB)\n', ...
    KK{iR1}.edge.SINR_p5, KK{iF3}.edge.SINR_p5, KK{iF3}.edge.SINR_p5-KK{iR1}.edge.SINR_p5);
fprintf('[H1] Penalizacion por interferencia (C/N - SINR) media: reuse1 = %.2f dB -> FFR(D=3) = %.2f dB\n', ...
    KK{iR1}.penalty_mean, KK{iF3}.penalty_mean);
fprintf('==================================================================================\n');

%% 5b. Desglose por ELEVACION del servidor (mecanismo que motiva la FFR ADAPTATIVA)
%  Con haces Earth-fixed y HPBW fijo, al bajar la elevacion el cluster de celdas se
%  COMPRIME angularmente visto desde el satelite -> los haces se solapan mas -> la
%  interferencia intra-satelite empeora. Es la dependencia que explotara ffr_policy
%  en H2 (mas proteccion de borde en geometria mala, mas eficiencia en geometria buena).
edgesEl = [25 35 45 60 90];
fprintf('\n--------- SINR_p5 (dB) POR TRAMO DE ELEVACION DEL SERVIDOR ---------\n');
fprintf('%-18s', 'esquema');
for e = 1:numel(edgesEl)-1, fprintf('%11s', sprintf('%d-%d deg',edgesEl(e),edgesEl(e+1))); end
fprintf('\n%s\n', repmat('-',1,18+11*(numel(edgesEl)-1)));
elBin = ctx.elev_deg;
for c = 1:nC
    fprintf('%-18s', cases{c}.name);
    for e = 1:numel(edgesEl)-1
        sel = ctx.cov & elBin >= edgesEl(e) & elBin < edgesEl(e+1);
        s   = FF{c}.SINR_dB(sel);
        if isempty(s), fprintf('%11s','-'); else, fprintf('%11.2f', prctile(s,5)); end
    end
    fprintf('\n');
end
fprintf('%-18s', 'muestras');
for e = 1:numel(edgesEl)-1
    fprintf('%11d', sum(ctx.cov(:) & elBin(:) >= edgesEl(e) & elBin(:) < edgesEl(e+1)));
end
% Celdas OCUPADAS por la rejilla: comprueba que la asociacion usuario-haz NO se
% degrada con la compresion angular. Geometricamente los usuarios (radio 40 km ~
% 3.2 celdas) deben caer en ~35-45 celdas de las 91; si a baja elevacion se
% dispersaran por las 91, la asociacion estaria rota. (No lo esta: se mantiene.)
fprintf('\n%-18s', 'celdas ocupadas');
for e = 1:numel(edgesEl)-1
    nOcc = [];
    for k = 1:ctx.Nt
        sel = ctx.cov(:,k) & ctx.elev_deg(:,k) >= edgesEl(e) & ctx.elev_deg(:,k) < edgesEl(e+1);
        if any(sel), nOcc(end+1) = numel(unique(ctx.beam(sel,k))); end   %#ok<SAGROW>
    end
    if isempty(nOcc), fprintf('%11s','-'); else, fprintf('%11.0f', mean(nOcc)); end
end
fprintf('   (de %d haces)\n%s\n', BL.nBeams, repmat('-',1,18+11*(numel(edgesEl)-1)));
% OJO: a densidad REAL los tramos de elevacion BAJA pueden quedar VACIOS (el mejor
% de 1584 satelites no baja de ~49 deg), asi que el diagnostico de compresion
% angular solo se imprime si hay muestras a ambos lados. Es en si mismo el
% resultado: el mecanismo existe (se mide en E3b, donde el eje SI se recorre) pero
% a densidad real el sistema no visita la region donde muerde.
loEl = ctx.cov & elBin < 35;
hiEl = ctx.cov & elBin >= 60;
if any(loEl(:)) && any(hiEl(:))
    fprintf(['[diagnostico] offset angular medio al haz propio: %.3f deg a <35 deg de elev. ' ...
             'vs %.3f deg a >60 deg (s=%.3f deg).\n'], ...
        mean(ctx.theta_deg(loEl)), mean(ctx.theta_deg(hiEl)), BL.spacing_deg);
    fprintf(['[MECANISMO] La asociacion usuario-haz se MANTIENE (celdas ocupadas ~ constante),\n' ...
             '            pero al comprimirse el cluster los haces vecinos se SOLAPAN mucho mas\n' ...
             '            (offset propio %.2f deg frente a un HPBW de %.2f deg) y el reuso-1 se hunde.\n' ...
             '            El plan de celdas Earth-fixed dimensionado a nadir solo es EFICIENTE a\n' ...
             '            elevacion alta: es la region de validez y la que motiva la FFR adaptativa.\n'], ...
        mean(ctx.theta_deg(loEl)), cfg.radio.beamwidth3dB_deg);
else
    fprintf(['[diagnostico] offset angular medio al haz propio a >60 deg de elevacion: %.3f deg ' ...
             '(s=%.3f deg, HPBW=%.2f deg).\n'], ...
        mean(ctx.theta_deg(hiEl)), BL.spacing_deg, cfg.radio.beamwidth3dB_deg);
    fprintf(['[MECANISMO] Con la constelacion REAL no hay muestras por debajo de 35 deg de\n' ...
             '            elevacion (elev minima observada %.2f deg): el mejor de %d satelites\n' ...
             '            mantiene al usuario en geometria favorable de forma permanente, luego\n' ...
             '            la COMPRESION ANGULAR del cluster Earth-fixed no llega a actuar. El\n' ...
             '            mecanismo se cuantifica en E3b, que baja la mascara a proposito para\n' ...
             '            recorrer el eje; aqui su AUSENCIA es el resultado.\n'], ...
        min(ctx.elev_deg(ctx.cov)), sats.N);
end

%% 5c. HEADLINE con GEOMETRIA VALIDA (elevacion >= elMin_valid)
%  Es la comparacion que se defiende: el plan de celdas Earth-fixed esta
%  dimensionado para vista cercana al nadir, luego los KPIs se evaluan donde el
%  cluster no esta comprimido. CONTRASTE CON E0: en el tramo 60-90 deg el reuso-1
%  da SINR_p5 ~ -4.2 dB frente a la CIR_p5 = -2.28 dB de la calibracion E0; la
%  diferencia es coherente con tener 91 haces co-canal en lugar de 19 y usuarios
%  repartidos en 3 anillos de celdas en lugar de solo el haz central.
elMin_valid = 45;
maskV = ctx.cov & ctx.elev_deg >= elMin_valid;
KV = cell(1,nC);
for c = 1:nC, KV{c} = compute_kpis(cfg, FF{c}, AA{c}, maskV); end

% ¿MUERDE la mascara? A densidad REAL (T=1584) el servidor es el mejor de 1584
% satelites y su elevacion nunca baja de ~52 deg (hallazgo de diseno de E3b), luego
% el filtro de "geometria valida" NO descarta ninguna muestra y las dos tablas
% coinciden: la de la ventana completa YA es la defendible. A T=66 si mordia
% (descartaba ~76% de las muestras) y por eso la tabla filtrada era la de cabecera.
% El bloque se conserva para poder repetir esa comparacion con constelacion reducida.
nCov  = sum(ctx.cov(:));
nVal  = sum(maskV(:));
elMinObs = min(ctx.elev_deg(ctx.cov));
maskBinds = nVal < nCov;
if maskBinds, bindTxt = 'MUERDE'; else, bindTxt = 'NO MUERDE (las dos tablas coinciden)'; end
fprintf(['\n[mascara] geometria valida (elev >= %d deg): %d de %d muestras con cobertura ' ...
         '(%.1f%%) | elevacion MINIMA observada del servidor = %.2f deg -> la mascara %s\n'], ...
    elMin_valid, nVal, nCov, 100*nVal/max(nCov,1), elMinObs, bindTxt);

fprintf('\n=========== TABLA CON GEOMETRIA VALIDA (elevacion >= %d deg, %d muestras) ===========\n', ...
    elMin_valid, sum(maskV(:)));
fprintf('%-18s %8s %9s %9s %9s %9s %7s %9s\n', ...
    'esquema','SINR_p5','R_p5','R_media','R_agreg','P(cob)','Jain','R_p5_BORDE');
fprintf('%-18s %8s %9s %9s %9s %9s %7s %9s\n', ...
    '','[dB]','[Mbps]','[Mbps]','[Gbps]','','','[Mbps]');
fprintf('%s\n', repmat('-',1,84));
for c = 1:nC
    K = KV{c};
    fprintf('%-18s %8.2f %9.1f %9.1f %9.2f %9.2f %7.3f %9.1f\n', ...
        cases{c}.name, K.all.SINR_p5, K.all.R_p5_Mbps, K.all.R_mean_Mbps, ...
        K.all.R_agg_Mbps/1e3, K.all.Pcov, K.all.jain, K.edge.R_p5_Mbps);
end
fprintf('%s\n', repmat('-',1,84));
fprintf('[H1 | geometria valida] R_p5 de BORDE: reuse1 = %.1f -> FFR(D=3) = %.1f Mbps (x%.2f)\n', ...
    KV{iR1}.edge.R_p5_Mbps, KV{iF3}.edge.R_p5_Mbps, ...
    KV{iF3}.edge.R_p5_Mbps/max(KV{iR1}.edge.R_p5_Mbps,eps));
fprintf('[H1 | geometria valida] SINR_p5 de BORDE: reuse1 = %.2f -> FFR(D=3) = %.2f dB (%+.2f dB)\n', ...
    KV{iR1}.edge.SINR_p5, KV{iF3}.edge.SINR_p5, KV{iF3}.edge.SINR_p5-KV{iR1}.edge.SINR_p5);

%% 5d. CRITERIO DE VIABILIDAD (H1) -- salida de DECISION
%  Formalizacion del umbral de viabilidad (antes tarea abierta). El criterio vive
%  en el MOTOR (cfg.viab + compute_kpis -> K.viab), no en este runner: E3 lo
%  reutilizara tal cual para barrer la densidad y localizar el MARGEN DE
%  VIABILIDAD = densidad maxima con verdict = 'viable'.
%
%  TRAZABILIDAD: percentil 5 de BORDE = "cell edge" de ITU-R M.2135-1 / M.2410-0,
%  heredado por 3GPP TR 38.821 (la misma metodologia con la que se calibro E0);
%  suelo gamma_floor = umbral de demodulacion del MCS mas robusto de 5G NR (QPSK,
%  10% BLER, TS 38.214). El criterio es INDEPENDIENTE DE LA BANDA: la SINR ya
%  incorpora lo que depende de la portadora (FSPL, lluvia P.618) y el umbral lo
%  fija la modulacion, luego se aplica en Ku sin reajuste aunque E0 se calibrara
%  con los parametros Ka de TR 38.821.
vb = cfg.viab;
fprintf('\n============================ CRITERIO DE VIABILIDAD (H1) ============================\n');
fprintf(['metrica = %s | suelo QPSK 5G NR = %.1f dB (TS 38.214) | umbral = %.1f dB | ' ...
         'cobertura exigida = %.0f%% | SE_p5 objetivo = %.2f b/s/Hz (M.2410-0 Rural)\n'], ...
    vb.metric, vb.gamma_floor_dB, vb.gamma_th_dB, 100*vb.cov_target, vb.se_target_bhz);
fprintf('poblacion juzgada: usuarios de BORDE (~ALLOC.isCenter), IDENTICA en todos los esquemas\n');

vtabs = {KK, sprintf('VENTANA COMPLETA (elev >= %g deg, %d muestras de borde)', ...
                     cfg.geom.minElev, KK{1}.viab.n_edge); ...
         KV, sprintf('GEOMETRIA VALIDA (elev >= %d deg, %d muestras de borde)', ...
                     elMin_valid, KV{1}.viab.n_edge)};
for t = 1:size(vtabs,1)
    KT = vtabs{t,1};
    fprintf('\n--- %s ---\n', vtabs{t,2});
    fprintf('%-18s %11s %11s %10s %11s %11s  %s\n', ...
        'esquema','SINR_p5_B','SE_p5_B','cobertura','margen_suelo','margen_umbr','VEREDICTO');
    fprintf('%-18s %11s %11s %10s %11s %11s\n', ...
        '','[dB]','[b/s/Hz]','P(>=umbr)','[dB]','[dB]');
    fprintf('%s\n', repmat('-',1,90));
    seTxt = {'  (SE_p5 objetivo: OK)', '  (SE_p5 objetivo: NO)'};
    for c = 1:nC
        V = KT{c}.viab;
        fprintf('%-18s %11.2f %11.3f %10.2f %11.2f %11.2f  %-9s%s\n', ...
            cases{c}.name, V.SINR_edge_p5, V.SE_edge_p5, V.coverage, ...
            V.margin_floor_dB, V.margin_th_dB, upper(V.verdict), ...
            seTxt{2 - double(V.meets_se_target)});
    end
    fprintf('%s\n', repmat('-',1,90));
end

fprintf(['\n[H1 | viabilidad] reuse1 = %s (SINR_p5 borde %.2f dB, cobertura %.2f) -> ' ...
         'FFR(D=3) = %s (%.2f dB, %.2f) | FFR(D=4) = %s (%.2f dB, %.2f)   [geometria valida]\n'], ...
    upper(KV{iR1}.viab.verdict), KV{iR1}.viab.SINR_edge_p5, KV{iR1}.viab.coverage, ...
    upper(KV{iF3}.viab.verdict), KV{iF3}.viab.SINR_edge_p5, KV{iF3}.viab.coverage, ...
    upper(KV{iF4}.viab.verdict), KV{iF4}.viab.SINR_edge_p5, KV{iF4}.viab.coverage);
fprintf(['[E3] Este veredicto es la SALIDA DE DECISION del barrido de saturacion: la\n' ...
         '     densidad maxima de satelites con verdict = ''viable'' es el MARGEN DE\n' ...
         '     VIABILIDAD, resultado central de H1.\n']);
fprintf('=====================================================================================\n');

%% 5e. DONDE CAE EL PERCENTIL 5 GLOBAL: en el CENTRO o en el BORDE?
%  HIPOTESIS a contrastar: en los esquemas FFR el percentil 5 GLOBAL lo fijan los
%  usuarios de CENTRO, no los de borde. El motivo seria estructural: el interior
%  de la FFR usa la sub-banda plena, REUTILIZADA POR LOS 91 HACES (es decir, la
%  SINR de reuso-1, sin ninguna proteccion) y ademas con solo Bc = alpha*B de
%  ancho, mientras el borde tiene reuso-Delta. Si se confirma, explica dos cosas
%  que en la TABLA COMPARATIVA parecen coincidencia:
%    (1) FFR(D=3) y FFR(D=4) dan el MISMO SINR_p5 y el MISMO R_p5 GLOBALES,
%        porque la clase que fija el percentil (centro) NO depende de Delta; y
%    (2) la FFR "pierde" en R_p5 global frente a reuso-1 pese a ganar en el borde.
%  Es decir: la FFR no elimina el cuello de botella, lo DESPLAZA del borde al
%  centro. Es material para la defensa porque el KPI reina (p5 de BORDE) y el p5
%  GLOBAL dejan de medir a la misma poblacion.
%
%  Se mide sobre la MISMA submuestra que compute_kpis (valid = SINR y B_user
%  finitos y B_user > 0), para que el percentil sea EXACTAMENTE el reportado.
p5split = repmat(struct('name','','thR_Mbps',NaN,'fracC_R',NaN,'liftC_R',NaN, ...
                        'thS_dB',NaN,'fracC_S',NaN,'nBelowR',0), 1, nC);
fprintf('\n============ DONDE CAE EL PERCENTIL 5 GLOBAL (centro vs borde) ============\n');
fprintf(['Clasificacion centro/borde IDENTICA en los %d esquemas (verificado en 4b), luego\n' ...
         'las fracciones son comparables entre esquemas. LINEA BASE: %.1f%% de las muestras\n' ...
         'validas son de CENTRO, asi que un reparto "neutro" daria ese mismo %.1f%% bajo el p5.\n'], ...
        nC, 100*KK{1}.centerFrac, 100*KK{1}.centerFrac);
fprintf('%-18s %9s %9s %9s | %9s %9s %9s\n', ...
    'esquema','p5 de R','%CENTRO','sobre-','p5 SINR','%CENTRO','sobre-');
fprintf('%-18s %9s %9s %9s | %9s %9s %9s\n', ...
    '','[Mbps]','bajo p5','rrepr.','[dB]','bajo p5','rrepr.');
fprintf('%s\n', repmat('-',1,74));
for c = 1:nC
    SINRc = FF{c}.SINR_dB;
    Buc   = AA{c}.B_user_Hz;
    isCc  = AA{c}.isCenter;
    validc = ~isnan(SINRc) & ~isnan(Buc) & (Buc > 0);

    Rc = nan(size(SINRc));
    Rc(validc) = Buc(validc) .* log2(1 + 10.^(SINRc(validc)/10)) / 1e6;   % Mbps

    baseC = sum(validc(:) & isCc(:)) / max(sum(validc(:)),1);   % fraccion de centro global

    % --- percentil 5 de THROUGHPUT
    thR    = prctile(Rc(validc), 5);
    belowR = validc & (Rc <= thR);
    fracCR = sum(belowR(:) & isCc(:)) / max(sum(belowR(:)),1);

    % --- percentil 5 de SINR
    thS    = prctile(SINRc(validc), 5);
    belowS = validc & (SINRc <= thS);
    fracCS = sum(belowS(:) & isCc(:)) / max(sum(belowS(:)),1);

    p5split(c) = struct('name',cases{c}.name, 'thR_Mbps',thR, 'fracC_R',fracCR, ...
                        'liftC_R',fracCR/max(baseC,eps), 'thS_dB',thS, ...
                        'fracC_S',fracCS, 'nBelowR',sum(belowR(:)));

    fprintf('%-18s %9.2f %8.1f%% %8.2fx | %9.2f %8.1f%% %8.2fx\n', ...
        cases{c}.name, thR, 100*fracCR, fracCR/max(baseC,eps), ...
        thS, 100*fracCS, fracCS/max(baseC,eps));
end
fprintf('%s\n', repmat('-',1,74));
fprintf(['LECTURA: "%%CENTRO bajo p5" = que fraccion de las muestras que caen POR DEBAJO del\n' ...
         'percentil 5 global son de CENTRO. "sobrerrepr." = esa fraccion dividida por la\n' ...
         'fraccion de centro en la poblacion (1.00x = neutro, 2.00x = el doble de lo que le\n' ...
         'tocaria). Si en los esquemas ffr sale ~100%% y ~2x, la hipotesis queda CONFIRMADA:\n' ...
         'el cuello de botella del p5 GLOBAL se ha desplazado del borde al CENTRO.\n']);
iF34 = [iF3 iF4];
fprintf(['[hipotesis] FFR(D=3) y FFR(D=4): SINR_p5 global = %.2f y %.2f dB (identicos si el p5\n' ...
         '            cae en el CENTRO, cuya SINR es la de reuso-1 y NO depende de Delta);\n' ...
         '            SINR_p5 de BORDE = %.2f y %.2f dB (ahi Delta SI actua).\n'], ...
    KK{iF34(1)}.all.SINR_p5, KK{iF34(2)}.all.SINR_p5, ...
    KK{iF34(1)}.edge.SINR_p5, KK{iF34(2)}.edge.SINR_p5);
fprintf('==========================================================================\n');

%% 6. (d) Barrido de alpha -> curva de compromiso
alphas = 0.1:0.1:0.9;
Dsweep = [3 4];
sw = struct('alpha',alphas,'Delta',Dsweep, ...
            'R_p5',nan(numel(Dsweep),numel(alphas)), ...
            'R_agg',nan(numel(Dsweep),numel(alphas)), ...
            'R_p5_edge',nan(numel(Dsweep),numel(alphas)), ...
            'R_agg_edge',nan(numel(Dsweep),numel(alphas)));

sw.R_p5_V      = nan(numel(Dsweep),numel(alphas));   % submuestra de geometria valida
sw.R_agg_V     = nan(numel(Dsweep),numel(alphas));
sw.R_p5_edge_V = nan(numel(Dsweep),numel(alphas));

fprintf('\n--------- BARRIDO DE ALPHA (curva de compromiso) ---------\n');
fprintf('%5s %6s | %10s %10s %12s | %10s %10s %12s\n', ...
    'Delta','alpha','R_p5','R_p5_BORDE','R_agreg','R_p5(V)','R_p5_B(V)','R_agreg(V)');
fprintf('%5s %6s | %10s %10s %12s | %10s %10s %12s\n', ...
    '','','[Mbps]','[Mbps]','[Gbps]','[Mbps]','[Mbps]','[Gbps]');
for id = 1:numel(Dsweep)
    for ia = 1:numel(alphas)
        cfg.ffr.scheme = 'ffr';
        cfg.ffr.Delta  = Dsweep(id);
        cfg.ffr.alpha  = alphas(ia);
        A = ffr_allocate(cfg, ctx, false);      % verbose = false: barrido silencioso
        F = compute_sinr_ffr(cfg, ctx, A);
        K = compute_kpis(cfg, F, A);
        Kv= compute_kpis(cfg, F, A, maskV);
        sw.R_p5(id,ia)       = K.all.R_p5_Mbps;
        sw.R_agg(id,ia)      = K.all.R_agg_Mbps;
        sw.R_p5_edge(id,ia)  = K.edge.R_p5_Mbps;
        sw.R_agg_edge(id,ia) = K.edge.R_agg_Mbps;
        sw.R_p5_V(id,ia)      = Kv.all.R_p5_Mbps;
        sw.R_agg_V(id,ia)     = Kv.all.R_agg_Mbps;
        sw.R_p5_edge_V(id,ia) = Kv.edge.R_p5_Mbps;
        fprintf('%5d %6.1f | %10.1f %10.1f %12.2f | %10.1f %10.1f %12.2f\n', ...
            Dsweep(id), alphas(ia), K.all.R_p5_Mbps, K.edge.R_p5_Mbps, K.all.R_agg_Mbps/1e3, ...
            Kv.all.R_p5_Mbps, Kv.edge.R_p5_Mbps, Kv.all.R_agg_Mbps/1e3);
    end
end

% LECTURA DEL COMPROMISO (honesta): R_agregado varia poco con alpha (la banda se
% conserva: B = Bc + Delta*Be), mientras R_p5 varia mucho. El optimo de compromiso
% lo fija por tanto R_p5, y el R_p5 de BORDE es MONOTONO DECRECIENTE en alpha
% (menos alpha = mas banda al borde = mas proteccion).
sw.alpha_best   = nan(1,numel(Dsweep));
sw.alpha_best_V = nan(1,numel(Dsweep));
for id = 1:numel(Dsweep)
    [~, ib]  = max(sw.R_p5(id,:));
    [~, ibv] = max(sw.R_p5_V(id,:));
    sw.alpha_best(id)   = alphas(ib);
    sw.alpha_best_V(id) = alphas(ibv);
    fprintf(['[H1/H2] Delta=%d: alpha de MEJOR COMPROMISO (max R_p5 global) = %.1f ' ...
             '(R_p5 = %.1f Mbps, R_agreg = %.2f Gbps; R_agreg varia solo %.1f%% en todo el barrido)\n'], ...
        Dsweep(id), alphas(ib), sw.R_p5(id,ib), sw.R_agg(id,ib)/1e3, ...
        100*(max(sw.R_agg(id,:))-min(sw.R_agg(id,:)))/mean(sw.R_agg(id,:)));
    fprintf(['        geometria valida (elev>=%d deg): alpha* = %.1f (R_p5 = %.1f Mbps) | ' ...
             'proteccion MAXIMA de borde en alpha=%.1f (R_p5_borde = %.1f Mbps)\n'], ...
        elMin_valid, alphas(ibv), sw.R_p5_V(id,ibv), alphas(1), sw.R_p5_edge_V(id,1));
end

%% 7. Figuras
%  Se exportan a PNG con save_fig (300 dpi, fuente unica) en FIGDIR. Los nombres
%  son DESCRIPTIVOS a proposito: en la memoria hay ~25 figuras y 'fig1/fig2' no
%  permite saber que muestra cada una sin abrirlas.

% (a) CDF de SINR por esquema
fA = figure('Name','CDF SINR por esquema','Color','w');
for c = 1:nC
    cdf = KK{c}.all.cdf_SINR;
    plot(cdf(:,1), cdf(:,2), 'LineWidth',1.6, 'DisplayName',cases{c}.name); hold on;
end
xline(cfg.kpi.gamma0_dB,'k--','DisplayName',sprintf('\\gamma_0 = %g dB',cfg.kpi.gamma0_dB));
grid on; xlabel('SINR [dB]'); ylabel('CDF');
title('CDF de SINR: reuso-1 vs reuso-\Delta vs FFR (layout multihaz Earth-fixed)');
legend('Location','southeast');
save_fig(fA, FIGDIR, 'e2_a_cdf_sinr_por_esquema');

% (b) CDF de throughput por esquema
fB = figure('Name','CDF throughput por esquema','Color','w');
for c = 1:nC
    cdf = KK{c}.all.cdf_R;
    semilogx(cdf(:,1), cdf(:,2), 'LineWidth',1.6, 'DisplayName',cases{c}.name); hold on;
end
grid on; xlabel('Throughput por usuario R [Mbps]'); ylabel('CDF');
title('CDF de throughput (Shannon): B_{user}\cdot log_2(1+SINR)');
legend('Location','southeast');
save_fig(fB, FIGDIR, 'e2_b_cdf_throughput_todos');

% (b2) CDF de throughput SOLO usuarios de BORDE (metrica reina de H1)
fC = figure('Name','CDF throughput BORDE','Color','w');
for c = 1:nC
    cdf = KK{c}.edge.cdf_R;
    if all(isnan(cdf(:))), continue; end
    semilogx(cdf(:,1), cdf(:,2), 'LineWidth',1.6, 'DisplayName',cases{c}.name); hold on;
end
grid on; xlabel('Throughput R [Mbps]'); ylabel('CDF');
title('CDF de throughput de los usuarios de BORDE (percentil 5 = metrica reina)');
legend('Location','southeast');
save_fig(fC, FIGDIR, 'e2_c_cdf_throughput_borde');

% (d) Curva de compromiso del barrido de alpha
fD = figure('Name','Compromiso alpha','Color','w');
yyaxis left
plot(alphas, sw.R_p5(1,:), 'o-','LineWidth',1.6,'DisplayName','R_{p5} (\Delta=3)'); hold on;
plot(alphas, sw.R_p5_edge(1,:), 's--','LineWidth',1.4,'DisplayName','R_{p5} borde (\Delta=3)');
ylabel('R_{p5} [Mbps]');
yyaxis right
plot(alphas, sw.R_agg(1,:)/1e3, '^-','LineWidth',1.6,'DisplayName','R_{agregado} (\Delta=3)');
ylabel('R_{agregado} [Gbps]');
grid on; xlabel('\alpha (fraccion de banda del interior)');
title('Curva de compromiso de la FFR: proteccion de borde vs eficiencia agregada');
legend('Location','best');
save_fig(fD, FIGDIR, 'e2_d_barrido_alpha_compromiso');

% (e) Mapa de celdas Earth-fixed, coloreado y rejilla de usuarios
fE = figure('Name','Celdas Earth-fixed y coloreado FFR','Color','w');
colD = ffr_coloring(BL, 3);
cmap = lines(3);
th = linspace(0,2*pi,40);
Rc = BL.spacing_km/sqrt(3);                       % circunradio de la celda hex
for b = 1:BL.nBeams
    p = BL.cell_km(b,:);
    fill(p(1)+Rc*cos(th), p(2)+Rc*sin(th), cmap(colD(b),:), ...
        'FaceAlpha',0.25, 'EdgeColor',[.6 .6 .6]); hold on;
end
plot(users.x_km, users.y_km, 'k.', 'MarkerSize',4);
axis equal; grid on; xlabel('Este [km]'); ylabel('Norte [km]');
title(sprintf('Celdas Earth-fixed (%d haces, s=%.1f km), coloreado \\Delta=3 y rejilla de usuarios', ...
    BL.nBeams, BL.spacing_km));
save_fig(fE, FIGDIR, 'e2_e_celdas_earthfixed_coloreado');

% (f) EL HALLAZGO DE 5e: donde cae el percentil 5 GLOBAL, centro vs borde.
%  Barras apiladas de la composicion de la cola (por debajo del p5 global) frente
%  a la LINEA BASE (la composicion de la poblacion). Si la barra de un esquema
%  coincide con la linea base, el p5 no discrimina entre clases; si es 100% de una
%  clase, el cuello de botella esta ENTERO en esa clase.
fF = figure('Name','Composicion del percentil 5 global','Color','w');
compC = 100*[p5split.fracC_R];              % % de CENTRO bajo el p5 de R
compB = 100 - compC;                        % % de BORDE
hb = bar([compC(:) compB(:)], 'stacked');
hb(1).FaceColor = [0.85 0.33 0.10];         % centro
hb(2).FaceColor = [0.00 0.45 0.74];         % borde
hold on;
base = 100*KK{1}.centerFrac;
yline(base, 'k--', 'LineWidth',1.6, ...
    'Label',sprintf('linea base: %.1f%% de centro en la poblacion', base), ...
    'LabelHorizontalAlignment','left');
set(gca,'XTickLabel',cellfun(@(c)c.name, cases, 'UniformOutput',false), ...
        'XTickLabelRotation',20);
ylim([0 100]); ylabel('composicion de las muestras BAJO el percentil 5 global [%]');
grid on;
legend({'usuarios de CENTRO','usuarios de BORDE','linea base (poblacion)'}, ...
       'Location','eastoutside');
title({'Donde cae el percentil 5 GLOBAL: la FFR DESPLAZA el cuello de botella al centro', ...
       'reuso-1 y reuso-\Delta: lo fija el BORDE   |   FFR: lo fija el CENTRO (100%)'});
save_fig(fF, FIGDIR, 'e2_f_composicion_del_percentil5_centro_vs_borde');

fprintf('\n[figuras] 6 PNG a 300 dpi en %s%s\n', FIGDIR, filesep);

%% 8. Guardar
% ctx.Grel / ctx.Glin son tensores [M x nBeams x Nt] de decenas de MB: se guarda
% una version LIGERA del contexto (todo lo demas se reconstruye reejecutando).
ctxLite = rmfield(ctx, {'Grel','Glin'});

timing = struct('warmup_p618_s', tWarm, 'experiment_s', toc(tExp), 'total_s', toc(tRun));
fprintf('\n[tiempos] experimento %.1f s (+ %.1f s de warmup P.618) | total %.1f s\n', ...
    timing.experiment_s, timing.warmup_p618_s, timing.total_s);

save('ffr_results.mat','cfg','BL','ctxLite','cases','KK','KV','AA','sw','alphas', ...
     'Dsweep','maskV','elMin_valid','p5split','timing','MEM','BINFO','maskBinds');
fprintf('\nResultados guardados en ffr_results.mat\n');
