function ALLOC = ffr_allocate(cfg, ctx, verbose)
%FFR_ALLOCATE  Clasificacion centro/borde y particion estricta de la banda (FFR).
%   ALLOC = ffr_allocate(cfg, ctx)
%   ALLOC = ffr_allocate(cfg, ctx, verbose)    (verbose = false: sin log; util en
%                                               barridos de alpha)
%
%   PIEZA DE MOTOR de la FFR. Decide, para cada usuario m e instante k:
%     (1) si es de CENTRO o de BORDE de su celda,
%     (2) que ANCHO DE BANDA se le asigna,
%     (3) en que SUB-BANDA transmite (lo que fija quien le interfiere).
%
%   PARTICION ESTRICTA (sin solape entre sub-bandas), con B = radio_band_Hz(cfg):
%       Bc = alpha*B                     interior (reutilizado por TODOS los haces)
%       Be = (1-alpha)*B/Delta           una sub-banda de borde por color
%       B  = Bc + Delta*Be               (identidad que se verifica)
%
%   CLASIFICADOR (cfg.ffr.classifier):
%     'sinr'  - CENTRO si la SINR de referencia en REUSO-1 >= tau_dB. Es el
%               criterio operativo (mide la calidad realmente sufrida, integra
%               geometria + interferencia). Si ctx.SINR_ref_dB esta vacio se
%               calcula aqui internamente con una asignacion reuso-1 (una sola vez;
%               conviene precalcularlo en el runner y pasarlo en ctx para reusarlo).
%     'theta' - CENTRO si el offset angular del usuario a su centro de haz
%               <= tau_theta_deg. Criterio puramente geometrico, mas simple de
%               justificar analiticamente y sin realimentacion del receptor.
%
%   ESQUEMAS (cfg.ffr.scheme). Toda la diferencia entre esquemas queda codificada
%   en ALLOC (la regla de co-canal de compute_sinr_ffr es UNICA: subband 0 => todos
%   los haces son co-canal; subband c>0 => solo los haces de color c):
%     'reuse1' - todos los usuarios usan B; sub-banda 0 -> co-canal = TODOS los haces.
%     'reuseD' - todos los usuarios usan B/Delta; sub-banda = color del haz
%                -> co-canal = solo los haces del MISMO color.
%     'ffr'    - CENTRO: Bc, sub-banda 0 (co-canal = todos, pero el centro esta
%                bien iluminado y lo tolera).
%                BORDE : Be, sub-banda = color del haz (co-canal = mismo color).
%
%   Salidas (struct ALLOC):
%     .isCenter   [M x Nt] logico (false tambien donde no hay cobertura)
%     .B_user_Hz  [M x Nt] ancho asignado [Hz]  (NaN sin cobertura)
%     .subband    [M x Nt] 0 = interior/banda plena ; 1..Delta = color de borde
%                          (NaN sin cobertura)
%     .color      [nBeams x 1] color de cada haz (de ffr_coloring)
%     .scheme, .Delta, .Bc_Hz, .Be_Hz, .B_Hz
%     .tau, .alpha  [1 x Nt] valores usados (de ffr_policy)
%     .POL          diagnostico de la politica (gancho H2)
%     .COL          diagnostico del coloreado (distancia minima entre co-canales)
%     .centerFrac   [1 x Nt] fraccion de usuarios de centro (entre los cubiertos)
%
%   Aportacion propia del TFG (logica de reuso fraccionario).

if nargin < 3 || isempty(verbose), verbose = true; end

M   = ctx.M;
Nt  = ctx.Nt;
cov = ctx.cov;
B   = radio_band_Hz(cfg);              % derivado de B_MHz en el punto de uso (G3)

scheme = lower(cfg.ffr.scheme);
Delta  = cfg.ffr.Delta;
if strcmp(scheme,'reuse1'), Delta = 1; end     % reuso-1 = un solo color

%% 1. SINR de REFERENCIA en reuso-1 (la usan el clasificador 'sinr' y, en modo
%    'quantile', tambien la politica: por eso se calcula ANTES de ffr_policy).
%    No se duplica fisica: se reutiliza compute_sinr_ffr con una asignacion de
%    reuso pleno. Conviene precalcularla en el runner y pasarla en ctx.SINR_ref_dB
%    para que TODOS los esquemas clasifiquen con la MISMA referencia.
if strcmpi(cfg.ffr.classifier,'sinr') && isempty(ctx.SINR_ref_dB)
    REF = struct('subband', zeros(M,Nt), 'B_user_Hz', B*ones(M,Nt), ...
                 'color', ones(ctx.nBeams,1), 'Delta', 1);
    REF.subband(~cov)   = NaN;
    REF.B_user_Hz(~cov) = NaN;
    tmp = compute_sinr_ffr(cfg, ctx, REF);
    ctx.SINR_ref_dB = tmp.SINR_dB;
end

%% 2. Politica: (tau, alpha) por instante  <-- gancho de la FFR adaptativa (H2)
[tau, alpha, POL] = ffr_policy(cfg, ctx);

%% 3. Coloreado del layout (que haces son co-canal entre si)
[color, COL] = ffr_coloring(ctx.BL, Delta, verbose);

%% 4. Clasificacion centro / borde
switch lower(cfg.ffr.classifier)
    case 'sinr'
        isCenter = cov & (ctx.SINR_ref_dB >= reshape(tau,1,Nt));
    case 'theta'
        isCenter = cov & (ctx.theta_deg <= reshape(tau,1,Nt));
    otherwise
        error('ffr_allocate:classifier', ...
            'Clasificador desconocido: %s (usa ''sinr'' o ''theta'').', cfg.ffr.classifier);
end

%% 5. Particion estricta de la banda (por instante, porque alpha puede variar)
Bc = alpha .* B;                       % [1 x Nt]  interior
Be = (1 - alpha) .* B ./ Delta;        % [1 x Nt]  cada sub-banda de borde
% Identidad B = Bc + Delta*Be (comprobacion de cordura del modulo)
resid = abs(Bc + Delta*Be - B);
if any(resid > 1e-6 * B)
    error('ffr_allocate:partition', ...
        'La particion no conserva la banda: max|Bc+Delta*Be-B| = %.3g Hz.', max(resid));
end

%% 6. Sub-banda de cada usuario segun el esquema (define QUIEN le interfiere)
colUser = nan(M, Nt);                              % color del haz de cada usuario
colUser(cov) = color(ctx.beam(cov));

subband = nan(M, Nt);
Bsub    = nan(M, Nt);                              % ancho de la SUB-BANDA usada
isEdge  = cov & ~isCenter;

switch scheme
    case 'reuse1'
        subband(cov) = 0;                          % banda plena: co-canal = todos
        Bsub(cov)    = B;
    case 'reused'
        subband(cov) = colUser(cov);               % co-canal = mismo color
        Bsub(cov)    = B / Delta;
    case 'ffr'
        BcM = repmat(Bc, M, 1);   BeM = repmat(Be, M, 1);
        subband(isCenter) = 0;                     % interior: co-canal = todos
        Bsub(isCenter)    = BcM(isCenter);
        subband(isEdge)   = colUser(isEdge);       % borde: co-canal = mismo color
        Bsub(isEdge)      = BeM(isEdge);
    otherwise
        error('ffr_allocate:scheme', ...
            'Esquema desconocido: %s (usa ''reuse1'', ''reuseD'' o ''ffr'').', cfg.ffr.scheme);
end

%% 7. CONVENCION DE PLANIFICACION -> ancho por USUARIO
%  Esta eleccion es METODOLOGICAMENTE CRITICA para que la comparacion sea limpia:
%
%  'peak'  - cada usuario recibe la sub-banda COMPLETA (tasa "asignable de pico").
%            SESGA la comparacion CONTRA la FFR: reuso-1 y reuso-Delta dan a todos
%            los usuarios de la celda la misma sub-banda, mientras la FFR PARTE la
%            poblacion en dos grupos y cada grupo solo puede usar su trozo, luego
%            la FFR paga el reparto sin que se le reconozca que su sub-banda la
%            comparten MENOS usuarios.
%  'share' - (POR DEFECTO) reparto equitativo de cada sub-banda entre los usuarios
%            de la celda que la usan (round-robin, buffers llenos), que es la
%            convencion habitual en los estudios de FFR:
%                B_user = B_sub / (n0 * f_clase)
%            con n0 = carga nominal de usuarios por celda (media de usuarios con
%            cobertura por celda ocupada, medida en la propia rejilla) y f_clase la
%            fraccion de usuarios de esa clase. Con esta convencion el caudal
%            AGREGADO por celda es Bc*log2(1+SINR_c) + Be*log2(1+SINR_b), es decir
%            NO depende de n0: la carga solo escala el eje de tasas por usuario,
%            identica en todos los esquemas, y la comparacion queda limpia.
%            Se usa la carga MEDIA (campo medio) y no el recuento celda a celda
%            para no arrastrar el artefacto de las celdas del borde de la rejilla,
%            que estan parcialmente pobladas.
sched = 'share';
if isfield(cfg.ffr,'sched') && ~isempty(cfg.ffr.sched), sched = lower(cfg.ffr.sched); end

B_user = nan(M, Nt);
switch sched
    case 'peak'
        B_user = Bsub;
        n0 = ones(1, Nt);
    case 'share'
        % Carga nominal por celda en cada instante: usuarios cubiertos / celdas ocupadas
        n0 = ones(1, Nt);
        for k = 1:Nt
            bk = ctx.beam(cov(:,k), k);
            if isempty(bk), continue; end
            n0(k) = numel(bk) / numel(unique(bk));
        end
        fC  = sum(isCenter,1) ./ max(sum(cov,1), 1);      % [1 x Nt] fraccion de centro
        fE  = 1 - fC;
        N0m = repmat(n0, M, 1);
        switch scheme
            case {'reuse1','reused'}
                B_user(cov) = Bsub(cov) ./ N0m(cov);
            case 'ffr'
                NCm = repmat(max(n0 .* fC, eps), M, 1);   % usuarios que comparten Bc
                NEm = repmat(max(n0 .* fE, eps), M, 1);   % usuarios que comparten Be
                B_user(isCenter) = Bsub(isCenter) ./ NCm(isCenter);
                B_user(isEdge)   = Bsub(isEdge)   ./ NEm(isEdge);
        end
    otherwise
        error('ffr_allocate:sched', ...
            'cfg.ffr.sched desconocido: %s (usa ''share'' o ''peak'').', sched);
end

%% 8. Salidas
ALLOC.isCenter  = isCenter;
ALLOC.SINR_ref_dB = ctx.SINR_ref_dB;   % referencia usada (para que el runner la cachee)
ALLOC.B_user_Hz = B_user;
ALLOC.B_sub_Hz  = Bsub;                % ancho de la sub-banda (antes del reparto)
ALLOC.sched     = sched;
ALLOC.load_n0   = n0;                  % carga nominal de usuarios por celda
ALLOC.subband   = subband;
ALLOC.color     = color;
ALLOC.scheme    = cfg.ffr.scheme;
ALLOC.Delta     = Delta;
ALLOC.B_Hz      = B;
ALLOC.Bc_Hz     = Bc;
ALLOC.Be_Hz     = Be;
ALLOC.tau       = tau;
ALLOC.alpha     = alpha;
ALLOC.POL       = POL;
ALLOC.COL       = COL;
ALLOC.centerFrac = sum(isCenter,1) ./ max(sum(cov,1), 1);

if verbose
    nCov = sum(cov(:));
    fprintf(['[ffr_allocate] esquema=%-6s Delta=%d alpha=%.2f tau=%.2f (%s) | ' ...
             'Bc=%.1f MHz Be=%.1f MHz | centro %.1f%% / borde %.1f%% | ' ...
             'sched=%s (n0=%.1f usuarios/celda)\n'], ...
        ALLOC.scheme, Delta, alpha(1), tau(1), lower(cfg.ffr.classifier), ...
        Bc(1)/1e6, Be(1)/1e6, ...
        100*sum(isCenter(:))/max(nCov,1), 100*(nCov-sum(isCenter(:)))/max(nCov,1), ...
        sched, mean(n0(sum(cov,1) > 0)));    % media SOLO sobre instantes con cobertura
end
end
