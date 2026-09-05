function [tau, alpha, POL] = ffr_policy(cfg, ctx)
%FFR_POLICY  Politica de la FFR: umbral de clasificacion y reparto de banda.
%   [tau, alpha]      = ffr_policy(cfg, ctx)
%   [tau, alpha, POL] = ffr_policy(cfg, ctx)
%
%   GANCHO DE LA FFR ADAPTATIVA (hipotesis H2). Es el UNICO punto del motor donde
%   se deciden los dos parametros libres de la FFR:
%
%     tau   - umbral del clasificador centro/borde. Sus unidades son las del
%             clasificador activo (cfg.ffr.classifier):
%               'sinr'  -> [dB]   (centro si SINR_ref reuso-1 >= tau)
%               'theta' -> [deg]  (centro si offset al centro de haz <= tau)
%     alpha - fraccion de la banda total reservada al INTERIOR (0..1):
%               B = alpha*B + Delta*((1-alpha)*B/Delta)
%
%   MODO ESTATICO (cfg.ffr.adaptive = false): devuelve las CONSTANTES de cfg,
%   replicadas por instante. La FFR es por tanto ESTATICA y sirve de linea base.
%
%   MODO ADAPTATIVO (cfg.ffr.adaptive = true, hipotesis H2 - EXPERIMENTO E5):
%   alpha varia POR INSTANTE con la ELEVACION del satelite servidor. La firma
%   devuelve [1 x Nt], asi que la adaptacion NO toca ffr_allocate,
%   compute_sinr_ffr, compute_kpis ni los runners.
%
%   REGLA (interpolacion lineal SATURADA entre dos anclas de elevacion):
%
%       w(t)     = clip( (elev(t) - elMin) / (elMax - elMin), 0, 1 )
%       alpha(t) = alpha_min + (alpha_max - alpha_min) * w(t)        [adapt_dir = +1]
%       alpha(t) = alpha_max + (alpha_min - alpha_max) * w(t)        [adapt_dir = -1]
%
%   con elev(t) = POL.elevMean_deg(t) (elevacion media del servidor entre los
%   usuarios con cobertura). alpha se decide POR INSTANTE y no por usuario porque
%   es una particion de banda a nivel de SISTEMA: todas las celdas comparten el
%   mismo plan Bc / Be en un instante dado.
%
%   SENTIDO DE LA ADAPTACION (adapt_dir) - punto delicado, leer antes de tocarlo:
%   en este motor **alpha es la fraccion del INTERIOR** (Bc = alpha*B), luego
%   MAS alpha = MENOS banda al BORDE. Por tanto "proteger el borde en geometria
%   mala" significa alpha BAJO a baja elevacion:
%
%     adapt_dir = +1 (DEFECTO, es el sentido MEDIDO en el experimento de la FFR
%       estatica): alpha CRECE con la elevacion.
%       - elevacion BAJA  -> alpha_min: el cluster Earth-fixed se comprime
%         angularmente, los haces se solapan y la sub-banda interior (reutilizada
%         por TODOS los haces) tiene SINR << 0 dB, luego darle banda no rinde:
%         se vuelca la banda en las sub-bandas de BORDE, que si estan protegidas.
%       - elevacion ALTA -> alpha_max: el interior tiene SINR buena, la sub-banda
%         interior rinde y se recupera eficiencia espectral agregada.
%       Evidencia: en el barrido de alpha fijo, alpha* = 0.5 en la ventana completa
%       frente a alpha* = 0.6 restringiendo a elevacion >= 45 deg (el optimo SUBE
%       con la elevacion).
%
%     adapt_dir = -1: sentido INVERSO (alpha DECRECE con la elevacion). Se deja
%       disponible para poder CONTRASTARLO empiricamente en E5 en lugar de
%       postularlo: si la adaptacion no fuese mas que un artefacto, ambos sentidos
%       rendirian igual. No es el sentido por defecto.
%
%   GANCHO PENDIENTE (extension futura, NO implementado): adaptar tambien el
%   umbral tau segun POL.edgeFrac_ref / POL.covFrac, es decir mover el CUANTIL de
%   corte centro/borde (cfg.ffr.tau_q) ademas del reparto de banda. Seria un
%   control en lazo cerrado sobre la CARGA de borde:
%
%       % q(t) = tau_q0 - kq * (POL.edgeFrac_ref(t) - edgeFrac_ref_nominal)*100;
%       % q(t) = min(max(q(t), cfg.ffr.tau_q_min), cfg.ffr.tau_q_max);
%       % -> recalcular tau0v(k) como el percentil (100-q(k)) de SINR_ref
%
%   Se deja SIN implementar a proposito: con tau por cuantil fijo (50%) las dos
%   clases tienen el mismo tamano en todos los esquemas y en todos los instantes,
%   de modo que la comparacion estatica/adaptativa aisla el efecto de ALPHA. Mover
%   tau a la vez cambiaria tambien la poblacion de borde sobre la que se mide
%   R_p5, y la ganancia ya no seria atribuible.
%
%   ctx (de ffr_context) debe traer: .elev_deg, .theta_deg, .cov, .Nt y, para el
%   clasificador 'sinr', .SINR_ref_dB.
%
%   Salidas:
%     tau, alpha - [1 x Nt]
%     POL        - struct de diagnostico del contexto y de la regla:
%       .adaptive        cfg.ffr.adaptive
%       .k               [1 x Nt] indice de instante
%       .elevMean_deg    [1 x Nt] elevacion media del servidor entre los usuarios
%                                 con cobertura (variable de control de la regla)
%       .edgeFrac_ref    [1 x Nt] fraccion de usuarios de BORDE con el umbral de
%                                 cfg (primera pasada; diagnostico de carga de borde)
%       .covFrac         [1 x Nt] fraccion de usuarios con cobertura
%       .classifier      clasificador activo
%       .tau_mode        'abs' | 'quantile'
%     Solo en modo adaptativo:
%       .w               [1 x Nt] peso saturado de la interpolacion (0..1)
%       .alpha_min/_max, .elMin_deg/.elMax_deg, .adapt_dir   parametros de la regla
%       .alpha_range     [min max] efectivamente recorrido por alpha(t)
%
%   Aportacion propia del TFG (politica de la FFR).

Nt = ctx.Nt;

%% 1. Diagnostico del contexto (lo que leera la regla adaptativa)
POL.adaptive   = cfg.ffr.adaptive;
POL.classifier = lower(cfg.ffr.classifier);
POL.k          = 1:Nt;
POL.covFrac    = sum(ctx.cov, 1) / max(ctx.M, 1);

elev = ctx.elev_deg;  elev(~ctx.cov) = NaN;
POL.elevMean_deg = mean(elev, 1, 'omitnan');

% Fraccion de borde de PRIMERA PASADA con el umbral de cfg (solo diagnostico:
% la clasificacion definitiva la hace ffr_allocate con el tau devuelto aqui).
switch POL.classifier
    case 'sinr'
        if isempty(ctx.SINR_ref_dB)
            POL.edgeFrac_ref = nan(1, Nt);
        else
            isEdge0 = ctx.cov & (ctx.SINR_ref_dB < cfg.ffr.tau_dB);
            POL.edgeFrac_ref = sum(isEdge0, 1) ./ max(sum(ctx.cov, 1), 1);
        end
    case 'theta'
        isEdge0 = ctx.cov & (ctx.theta_deg > cfg.ffr.tau_theta_deg);
        POL.edgeFrac_ref = sum(isEdge0, 1) ./ max(sum(ctx.cov, 1), 1);
    otherwise
        error('ffr_policy:classifier', ...
            'Clasificador desconocido: %s (usa ''sinr'' o ''theta'').', cfg.ffr.classifier);
end

%% 2. Decision de (tau, alpha)
switch POL.classifier
    case 'sinr',  tau0 = cfg.ffr.tau_dB;
    case 'theta', tau0 = cfg.ffr.tau_theta_deg;
end
alpha0 = cfg.ffr.alpha;

tau_mode = 'abs';
if isfield(cfg.ffr,'tau_mode') && ~isempty(cfg.ffr.tau_mode), tau_mode = lower(cfg.ffr.tau_mode); end
POL.tau_mode = tau_mode;

if strcmp(tau_mode,'quantile')
    % --- tau POR INSTANTE para que la fraccion de CENTRO sea cfg.ffr.tau_q % ---
    % Motivo (importante): el nivel ABSOLUTO de SINR en reuso-1 depende mucho de la
    % geometria (a baja elevacion el cluster Earth-fixed se comprime angularmente y
    % la SINR se hunde), luego un tau fijo deja el 100% de usuarios a un lado y la
    % FFR degenera (todos borde -> peor que reuso-Delta; todos centro -> reuso-1).
    % Fijar el REPARTO en vez del nivel es tambien lo que hara la politica
    % adaptativa (H2): es un control en lazo cerrado sobre la carga de borde.
    q = cfg.ffr.tau_q;
    tau0v = nan(1, Nt);
    switch POL.classifier
        case 'sinr'
            % centro si SINR >= tau  ->  q% de centro => tau = percentil (100-q)
            X = ctx.SINR_ref_dB;
            if isempty(X)
                error('ffr_policy:noRef', ...
                    ['tau_mode = ''quantile'' con clasificador ''sinr'' exige ' ...
                     'ctx.SINR_ref_dB (SINR de referencia en reuso-1).']);
            end
            pq = 100 - q;
        case 'theta'
            % centro si theta <= tau  ->  q% de centro => tau = percentil q
            X  = ctx.theta_deg;
            pq = q;
    end
    for k = 1:Nt
        xk = X(ctx.cov(:,k), k);
        if ~isempty(xk), tau0v(k) = prctile(xk, pq); end
    end
    tau0v(isnan(tau0v)) = tau0;          % instantes sin cobertura: valor de cfg
    tau0 = tau0v;
end

tau = tau0 .* ones(1, Nt);      % tau: por ahora identico en ambos modos
                                % (el gancho adaptativo de tau queda pendiente,
                                %  ver la cabecera de la funcion)

if ~cfg.ffr.adaptive
    % --- FFR ESTATICA: alpha constante de configuracion ---
    alpha = alpha0 * ones(1, Nt);
else
    % --- FFR ADAPTATIVA (H2, experimento E5): alpha(t) segun la ELEVACION ---
    a_min  = getf(cfg.ffr, 'alpha_min', alpha0);
    a_max  = getf(cfg.ffr, 'alpha_max', alpha0);
    elMinP = getf(cfg.ffr, 'elMin_deg', cfg.geom.minElev);
    elMaxP = getf(cfg.ffr, 'elMax_deg', 90);
    adir   = getf(cfg.ffr, 'adapt_dir', +1);

    if a_min > a_max
        error('ffr_policy:alphaRange', ...
            'cfg.ffr.alpha_min (%.3f) debe ser <= cfg.ffr.alpha_max (%.3f).', a_min, a_max);
    end
    if elMaxP <= elMinP
        error('ffr_policy:elRange', ...
            'cfg.ffr.elMax_deg (%.2f) debe ser > cfg.ffr.elMin_deg (%.2f).', elMaxP, elMinP);
    end
    if ~ismember(adir, [-1 1])
        error('ffr_policy:adaptDir', 'cfg.ffr.adapt_dir debe ser +1 o -1 (recibido %g).', adir);
    end

    % Variable de control: elevacion media del servidor en cada instante.
    % Instantes SIN cobertura (elevMean = NaN) no tienen usuarios que asignar;
    % se les da alpha0 para dejar el vector bien definido (no afecta a ningun KPI).
    el = POL.elevMean_deg;
    w  = (el - elMinP) ./ (elMaxP - elMinP);
    w  = min(max(w, 0), 1);                       % saturacion en las anclas

    if adir > 0
        alpha = a_min + (a_max - a_min) .* w;     % alpha CRECE con la elevacion
    else
        alpha = a_max + (a_min - a_max) .* w;     % alpha DECRECE con la elevacion
    end

    bad = isnan(alpha);
    alpha(bad) = alpha0;
    w(bad)     = NaN;

    POL.w          = w;
    POL.alpha_min  = a_min;
    POL.alpha_max  = a_max;
    POL.elMin_deg  = elMinP;
    POL.elMax_deg  = elMaxP;
    POL.adapt_dir  = adir;
    POL.alpha_range = [min(alpha(~bad)), max(alpha(~bad))];

    % Cordura: si alpha_min == alpha_max la regla DEBE degenerar en la estatica.
    % Es la prueba de que la adaptacion no introduce ninguna ventaja espuria.
    if a_min == a_max && any(abs(alpha - a_min) > 0 & ~bad)
        error('ffr_policy:degenerate', ...
            'alpha_min == alpha_max pero alpha(t) no es constante: la regla esta mal.');
    end
end

% Comprobacion de cordura del reparto de banda
if any(alpha < 0 | alpha > 1)
    error('ffr_policy:alpha', 'alpha debe estar en [0,1] (recibido min=%.3f max=%.3f).', ...
        min(alpha), max(alpha));
end
end

% -------------------------------------------------------------------------
function v = getf(s, name, defaultVal)
%GETF  Lectura de un campo opcional de cfg con valor por defecto (retrocompatible:
%      una cfg antigua sin los campos de E5 sigue funcionando).
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultVal;
end
end
