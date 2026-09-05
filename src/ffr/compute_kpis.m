function K = compute_kpis(cfg, FFR, ALLOC, mask)
%COMPUTE_KPIS  KPIs del sistema (capacidad Shannon) sobre usuarios x instantes.
%   K = compute_kpis(cfg, FFR, ALLOC)
%   K = compute_kpis(cfg, FFR, ALLOC, mask)   mask [M x Nt] logico: restringe la
%       submuestra (p.ej. solo instantes de geometria valida, elevacion >= X)
%
%   PIEZA DE MOTOR: convierte la SINR y el ancho asignado en las metricas con las
%   que se juzgan H1 y H2. Reutilizable por cualquier experimento (E0/E3/...) sin
%   cambios: solo depende de FFR (compute_sinr_ffr) y ALLOC (ffr_allocate).
%
%   CAPACIDAD (formula de Shannon, decision del tutor; tablas MODCOD 5G despues):
%       R(m,k) = B_user(m,k) * log2(1 + SINR_lin(m,k))          [bit/s]
%
%   Interpretacion de R (importante para la defensa): B_user es el ancho que
%   ffr_allocate asigna al usuario segun la convencion de planificacion
%   cfg.ffr.sched:
%     'share' (por defecto) - la sub-banda se reparte entre los usuarios de la celda
%              que la usan, luego R es el caudal REAL por usuario y R_agregado es el
%              caudal agregado de la poblacion. Es la convencion JUSTA para comparar
%              esquemas: la FFR parte la poblacion y su sub-banda la comparten menos
%              usuarios, algo que la convencion 'peak' no le reconoce.
%     'peak'  - R es la tasa asignable de PICO (sub-banda completa). Util como
%              diagnostico, pero sesga la comparacion contra la FFR.
%   En cualquiera de las dos la comparacion es legitima porque la poblacion de
%   usuarios e instantes es EXACTAMENTE la misma en todos los esquemas.
%
%   COBERTURA: solo se promedia donde hay servicio (SINR no NaN y B_user > 0).
%   Ademas se reporta la fraccion de muestras con servicio (K.covFrac) para que
%   ningun esquema gane "por no servir".
%
%   Entradas:
%     cfg  - configuracion (usa cfg.kpi.gamma0_dB)
%     FFR  - salida de compute_sinr_ffr (usa FFR.SINR_dB, FFR.CN_dB, FFR.penalty_dB)
%     ALLOC- salida de ffr_allocate (usa ALLOC.B_user_Hz, ALLOC.isCenter)
%
%   Salidas (struct K):
%     .all, .center, .edge  - sub-structs con las mismas metricas, para el total y
%                             para el desglose CENTRO vs BORDE:
%         .n                numero de muestras (usuario,instante)
%         .SINR_mean/_p5/_p50/_p95        [dB]
%         .R_mean_Mbps/_p5_Mbps/_p50_Mbps/_p95_Mbps
%         .R_agg_Mbps       caudal agregado: media temporal de la suma sobre usuarios
%         .SE_mean_bpsHz    eficiencia espectral media log2(1+SINR)
%         .Pcov             P(SINR >= cfg.kpi.gamma0_dB)
%         .jain             indice de equidad de Jain sobre R
%         .cdf_SINR/.cdf_R  vectores [x, F] para las CDF
%     .covFrac              fraccion de muestras (usuario,instante) con servicio
%     .gamma0_dB            umbral de cobertura usado
%     .penalty_mean/_p5     penalizacion por interferencia C/N - SINR [dB] (H1)
%     .B_user_MHz_mean      ancho medio asignado [MHz]
%     .viab                 CRITERIO DE VIABILIDAD (H1), ver seccion siguiente
%
%   ---------------------------------------------------------------------------
%   CRITERIO DE VIABILIDAD (K.viab) - salida de DECISION de la hipotesis H1
%   ---------------------------------------------------------------------------
%   Formaliza el umbral de "viabilidad" del TFG. Parametros en cfg.viab.
%
%   POBLACION EVALUADA: los usuarios de BORDE, con EXACTAMENTE la misma
%   clasificacion que usa la FFR (ALLOC.isCenter -> borde = ~isCenter), es decir
%   la que ffr_allocate calcula ANTES de ramificar por esquema, con la SINR de
%   REFERENCIA en reuso-1 y el umbral tau de ffr_policy. Por construccion es
%   INDEPENDIENTE de cfg.ffr.scheme, de modo que reuse1, reuseD y ffr se juzgan
%   sobre la MISMA poblacion de borde y la comparacion es justa (reuse1 y reuseD
%   no tienen particion centro/borde nativa: se les aplica la de referencia).
%   La mascara es literalmente la misma que alimenta K.edge, y el modulo ABORTA
%   si K.viab.SINR_edge_p5 difiere de K.edge.SINR_p5 (comprobacion de cordura:
%   garantiza que el percentil NO se esta tomando sobre todos los usuarios).
%
%   TRAZABILIDAD (ITU-R -> 3GPP -> E0):
%     - Percentil 5 de borde: definicion normalizada de "cell edge" de
%       ITU-R M.2135-1 (5% user spectral efficiency, IMT-Advanced) y de
%       ITU-R M.2410-0 (5th percentile user spectral efficiency, IMT-2020),
%       heredada por 3GPP TR 38.821 en la metodologia de evaluacion NTN, que es
%       la MISMA con la que se calibro E0 en este TFG.
%     - Suelo de SINR (cfg.viab.gamma_floor_dB): umbral de demodulacion del MCS
%       mas robusto de 5G NR (QPSK con la tasa de codigo mas baja, criterio de
%       10% BLER; tablas MCS/CQI de 3GPP TS 38.214). Por debajo NO existe MODCOD
%       que cierre el enlace, aunque Shannon siga dando capacidad > 0.
%     - Objetivo de eficiencia espectral de borde (cfg.viab.se_target_bhz):
%       ITU-R M.2410-0 Tabla 1, escenario Rural eMBB (5th percentile user SE, DL).
%
%   INDEPENDENCIA DE BANDA (por que vale en Ku aunque E0 se calibrara en Ka):
%     el criterio se formula sobre la SINR de borde y sobre umbrales de
%     MODULACION, ambos independientes de la portadora. Lo que si depende de la
%     banda (FSPL, lluvia y gases ITU-R P.618, ganancias de apertura) ya esta
%     incorporado DENTRO de la SINR; una vez calculada, el umbral que decide si
%     el enlace cierra lo fija la modulacion/codificacion, no la frecuencia.
%     E0 valida el MOTOR (como se calcula la SINR), no la banda.
%
%   Campos de K.viab:
%     .metric            metrica usada (cfg.viab.metric; hoy 'sinr_edge_p5')
%     .SINR_edge_p5      [dB]        percentil 5 de la SINR de los usuarios de BORDE
%     .SE_edge_p5        [bit/s/Hz]  percentil 5 de log2(1+SINR_lin) de BORDE
%     .coverage          [0..1]      P(SINR_borde >= cfg.viab.gamma_th_dB)
%     .margin_floor_dB   [dB]        SINR_edge_p5 - gamma_floor_dB (margen al suelo
%                                    fisico; < 0 => ningun MODCOD cierra en el borde)
%     .margin_th_dB      [dB]        SINR_edge_p5 - gamma_th_dB (margen al umbral)
%     .verdict           'inviable' | 'marginal' | 'viable' | 'sin-datos'
%           'inviable' si SINR_edge_p5 <  gamma_floor_dB
%           'marginal' si gamma_floor_dB <= SINR_edge_p5 < gamma_th_dB
%           'viable'   si SINR_edge_p5 >= gamma_th_dB  Y  coverage >= cov_target
%           (si SINR_edge_p5 >= gamma_th_dB pero la cobertura no llega al objetivo,
%            el veredicto NO asciende a 'viable': queda en 'marginal')
%     .meets_se_target   logico, SE_edge_p5 >= cfg.viab.se_target_bhz. Se reporta
%                        aparte y NO entra en el veredicto: el veredicto se define
%                        sobre la SINR (magnitud del enlace), mientras que la SE de
%                        borde de M.2410-0 es un objetivo de servicio con el que se
%                        contrasta el resultado.
%     .n_edge            numero de muestras de borde evaluadas
%     .edge_from         trazabilidad de la clasificacion usada ('ALLOC.isCenter')
%     .gamma_floor_dB/.gamma_th_dB/.cov_target/.se_target_bhz  parametros aplicados

SINR = FFR.SINR_dB;
Bu   = ALLOC.B_user_Hz;

valid = ~isnan(SINR) & ~isnan(Bu) & (Bu > 0);
if nargin >= 4 && ~isempty(mask), valid = valid & mask; end

%% Capacidad Shannon
R_bps = nan(size(SINR));
R_bps(valid) = Bu(valid) .* log2(1 + 10.^(SINR(valid)/10));

%% Desglose
isC = ALLOC.isCenter & valid;
isE = ~ALLOC.isCenter & valid;

K.all    = kpi_block(SINR, R_bps, valid, cfg.kpi.gamma0_dB);
K.center = kpi_block(SINR, R_bps, isC,   cfg.kpi.gamma0_dB);
K.edge   = kpi_block(SINR, R_bps, isE,   cfg.kpi.gamma0_dB);

%% Metricas globales
% Fraccion con servicio, referida a la submuestra evaluada (para que ningun
% esquema "gane" simplemente por dejar usuarios sin servir).
if nargin >= 4 && ~isempty(mask)
    K.covFrac = sum(valid(:)) / max(sum(mask(:)), 1);
else
    K.covFrac = sum(valid(:)) / numel(valid);
end
K.gamma0_dB = cfg.kpi.gamma0_dB;
K.scheme    = ALLOC.scheme;
K.Delta     = ALLOC.Delta;
K.alpha     = ALLOC.alpha(1);
K.centerFrac = sum(isC(:)) / max(sum(valid(:)), 1);

% Penalizacion por interferencia (metrica de H1: el "hueco" que abre el reuso)
pen = FFR.penalty_dB;
K.penalty_mean = mean(pen(valid));
K.penalty_p5   = prctile(pen(valid), 5);
K.penalty_p95  = prctile(pen(valid), 95);

K.B_user_MHz_mean = mean(Bu(valid)) / 1e6;

%% CRITERIO DE VIABILIDAD (H1) -- ver cabecera para fuentes y trazabilidad
% Se evalua sobre isE, la MISMA mascara de BORDE que alimenta K.edge (borde =
% ~ALLOC.isCenter, clasificacion de REFERENCIA comun a los tres esquemas).
K.viab = viability_block(cfg, SINR, isE);

% Comprobacion de cordura: el percentil de viabilidad DEBE ser el del bloque de
% borde. Si alguien cambiara la mascara (p.ej. tomando todos los usuarios), esto
% aborta en vez de devolver un veredicto silenciosamente equivocado.
if K.edge.n > 0 && abs(K.viab.SINR_edge_p5 - K.edge.SINR_p5) > 1e-9
    error('compute_kpis:viabMask', ...
        ['El bloque de viabilidad no esta evaluando la poblacion de BORDE ' ...
         '(SINR_edge_p5 = %.6f dB vs K.edge.SINR_p5 = %.6f dB).'], ...
        K.viab.SINR_edge_p5, K.edge.SINR_p5);
end
end

% -------------------------------------------------------------------------
function v = viability_block(cfg, SINR, isEdge)
%VIABILITY_BLOCK  Criterio de viabilidad de H1 sobre los usuarios de BORDE.
%   Fuentes: ITU-R M.2135-1 / M.2410-0 (percentil 5 de borde), 3GPP TS 38.214
%   (umbral de demodulacion QPSK, 10% BLER), 3GPP TR 38.821 (metodologia NTN,
%   la misma de la calibracion E0). Ver la cabecera de compute_kpis.

% Lectura defensiva de cfg.viab: una cfg antigua (p.ej. recuperada de un .mat
% anterior a este criterio) sigue funcionando con los valores por defecto.
V = struct('metric','sinr_edge_p5','gamma_floor_dB',-6.7,'gamma_th_dB',0.0, ...
           'cov_target',0.95,'se_target_bhz',0.12);
if isfield(cfg,'viab')
    f = fieldnames(V);
    for i = 1:numel(f)
        if isfield(cfg.viab, f{i}) && ~isempty(cfg.viab.(f{i}))
            V.(f{i}) = cfg.viab.(f{i});
        end
    end
end

% Hoy solo esta definida la metrica reina. El campo se valida (en vez de
% ignorarse) para que anadir otra metrica en el futuro sea un cambio explicito.
if ~strcmpi(V.metric, 'sinr_edge_p5')
    error('compute_kpis:viabMetric', ...
        'cfg.viab.metric = ''%s'' no soportada (hoy solo ''sinr_edge_p5'').', V.metric);
end

v.metric     = V.metric;
v.edge_from  = 'ALLOC.isCenter';   % trazabilidad: clasificacion de referencia de la FFR
v.n_edge     = sum(isEdge(:));

v.gamma_floor_dB = V.gamma_floor_dB;
v.gamma_th_dB    = V.gamma_th_dB;
v.cov_target     = V.cov_target;
v.se_target_bhz  = V.se_target_bhz;

if v.n_edge == 0
    v.SINR_edge_p5 = NaN;  v.SINR_edge_p1 = NaN;  v.SE_edge_p5 = NaN;  v.coverage = NaN;
    v.margin_floor_dB = NaN;  v.margin_th_dB = NaN;
    v.meets_se_target = false;
    v.verdict = 'sin-datos';
    return;
end

s = SINR(isEdge);                          % [dB] SOLO usuarios de BORDE
se = log2(1 + 10.^(s/10));                 % [bit/s/Hz] eficiencia espectral Shannon

v.SINR_edge_p5 = prctile(s, 5);            % <-- METRICA REINA del criterio
v.SINR_edge_p1 = prctile(s, 1);            % COLA (E4: eventos in-line inter-constelacion)
v.SE_edge_p5   = prctile(se, 5);
v.coverage     = mean(s >= V.gamma_th_dB);

v.margin_floor_dB = v.SINR_edge_p5 - V.gamma_floor_dB;
v.margin_th_dB    = v.SINR_edge_p5 - V.gamma_th_dB;
v.meets_se_target = v.SE_edge_p5 >= V.se_target_bhz;

% Veredicto de 3 niveles
if v.SINR_edge_p5 < V.gamma_floor_dB
    v.verdict = 'inviable';                % ningun MODCOD de 5G NR cierra en el borde
elseif v.SINR_edge_p5 >= V.gamma_th_dB && v.coverage >= V.cov_target
    v.verdict = 'viable';
else
    v.verdict = 'marginal';                % cierra con el MODCOD mas robusto, pero
end                                        % no alcanza umbral+cobertura exigidos
end

% -------------------------------------------------------------------------
function b = kpi_block(SINR, R_bps, mask, gamma0_dB)
%KPI_BLOCK  Estadisticos de un subconjunto de muestras (total / centro / borde).
b.n = sum(mask(:));
if b.n == 0
    f = {'SINR_mean','SINR_p1','SINR_p5','SINR_p50','SINR_p95','R_mean_Mbps', ...
         'R_p1_Mbps','R_p5_Mbps','R_p50_Mbps','R_p95_Mbps','R_agg_Mbps', ...
         'SE_mean_bpsHz','Pcov','jain'};
    for i = 1:numel(f), b.(f{i}) = NaN; end
    b.cdf_SINR = [NaN NaN];  b.cdf_R = [NaN NaN];
    return;
end

s = SINR(mask);
r = R_bps(mask) / 1e6;                       % Mbps

b.SINR_mean = mean(s);
b.SINR_p1   = prctile(s, 1);                  % COLA: la fija E4 (eventos in-line)
b.SINR_p5   = prctile(s, 5);
b.SINR_p50  = prctile(s, 50);
b.SINR_p95  = prctile(s, 95);

b.R_mean_Mbps = mean(r);
b.R_p1_Mbps   = prctile(r, 1);                % COLA
b.R_p5_Mbps   = prctile(r, 5);                % <-- METRICA REINA (throughput de borde)
b.R_p50_Mbps  = prctile(r, 50);
b.R_p95_Mbps  = prctile(r, 95);

% Caudal agregado: suma sobre usuarios en cada instante, promediada en el tiempo.
% (asi no depende del numero de instantes simulados)
Rk = R_bps;  Rk(~mask) = 0;
sumPerInstant = sum(Rk, 1) / 1e6;             % [1 x Nt] Mbps
nActive = sum(mask, 1) > 0;
b.R_agg_Mbps = mean(sumPerInstant(nActive));

b.SE_mean_bpsHz = mean(log2(1 + 10.^(s/10)));
b.Pcov          = mean(s >= gamma0_dB);

% Equidad de Jain sobre la tasa:  (sum R)^2 / (n * sum R^2)  in (0,1]
b.jain = sum(r)^2 / (numel(r) * sum(r.^2));

% CDFs para figuras
[Fs, Xs] = ecdf(s);   b.cdf_SINR = [Xs, Fs];
[Fr, Xr] = ecdf(r);   b.cdf_R    = [Xr, Fr];
end
