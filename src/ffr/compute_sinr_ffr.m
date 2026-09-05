function FFR = compute_sinr_ffr(cfg, ctx, ALLOC)
%COMPUTE_SINR_FFR  SINR por esquema de reuso sobre el layout multihaz Earth-fixed.
%   FFR = compute_sinr_ffr(cfg, ctx, ALLOC)
%
%   PIEZA DE MOTOR: nucleo de SINR de la FFR. Evalua la MISMA poblacion de usuarios
%   (ctx, de ffr_context) bajo la asignacion ALLOC (de ffr_allocate), de modo que
%   reuso-1, reuso-Delta y FFR son estrictamente comparables.
%
%   REGLA UNICA de conjunto co-canal (el esquema esta codificado en ALLOC.subband,
%   asi que aqui NO hay ningun switch por esquema):
%       ALLOC.subband(m,k) == 0  ->  co-canal = TODOS los haces del satelite salvo
%                                    el propio (banda plena / sub-banda interior,
%                                    que todos los haces reutilizan).
%       ALLOC.subband(m,k) == c  ->  co-canal = solo los haces de color c, salvo el
%                                    propio. Es la reduccion que aporta la FFR.
%
%   DESGLOSE INTRA / INTER SATELITE (se mantiene el del modulo existente):
%     - INTRA (dominante): haces co-canal del MISMO satelite servidor. El terminal
%       apunta al satelite, luego ve todos esos haces EN EJE (Grx_max comun a C e
%       I) y NO discrimina; FSPL y atmosfera tambien son comunes y se cancelan, asi
%       que I/C intra es un puro cociente de PATRONES (metodologia 3GPP TR 38.821).
%     - INTER (despreciable, medido en E0 y en run_interference_demo): se toma la
%       agregada de compute_interference (reuso-1, banda total) como DENSIDAD y se
%       reescala al ancho asignado y al reuso del esquema:
%           I_inter = Ipsd_inter + 10*log10(B_u[MHz]) + 10*log10(f)
%       con f = 1 en la sub-banda interior (todos los satelites la reutilizan) y
%       f = 1/Delta en una sub-banda de borde (un satelite interferente coincide en
%       color con probabilidad 1/Delta).  APROXIMACION documentada: es inmaterial
%       porque I_inter << N (~58 dB por debajo en E0), pero deja el desglose completo.
%
%   RUIDO: escala con el ANCHO REALMENTE ASIGNADO (imprescindible; fijarlo a la
%   banda total falsearia la comparacion entre esquemas):
%       N[dBW] = -228.6 + 10*log10(Tsys) + 10*log10(B_user[Hz])
%   Con la convencion de PSD constante de ffr_context, C tambien escala con B_user,
%   luego C/N es invariante al ancho y toda la ganancia de la FFR aparece via la
%   interferencia; el ancho entra despues en la capacidad (compute_kpis).
%
%   Salidas (struct FFR, [M x Nt], NaN sin cobertura):
%     .SINR_dB       SINR total (N + I_intra + I_inter)
%     .CN_dB         C/N del enlace servidor (sin interferencia)
%     .CIR_intra_dB  C/I intra-satelite (puro cociente de patrones)
%     .C_dBW, .N_dBW, .I_intra_dBW, .I_inter_dBW
%     .nCo           numero de haces co-canal (excluido el propio)
%     .penalty_dB    PENALIZACION por interferencia = C/N - SINR (metrica de H1)
%
%   Aportacion propia del TFG.

M   = ctx.M;
Nt  = ctx.Nt;
nB  = ctx.nBeams;
cov = ctx.cov;

Delta = ALLOC.Delta;
color = ALLOC.color;
Bu    = ALLOC.B_user_Hz;                    % [M x Nt] Hz
sb    = ALLOC.subband;                      % [M x Nt] 0..Delta

%% 1. Sumas de ganancia co-canal precalculadas (vectorizado, sin bucles)
%    sumAll(m,k)     = suma en lineal de la ganancia de TODOS los haces al usuario
%    sumCol(m,c,k)   = idem restringida a los haces de color c
Glin   = ctx.Glin;                                    % [M x nB x Nt]
sumAll = sum(Glin, 2, 'omitnan');                     % [M x 1 x Nt]
sumCol = zeros(M, Delta, Nt);
for c = 1:Delta
    sumCol(:,c,:) = sum(Glin(:, color == c, :), 2, 'omitnan');
end

%% 2. Ganancia del haz PROPIO y suma co-canal de cada usuario/instante
[mIdx, kIdx] = ndgrid(1:M, 1:Nt);
bOwn = ctx.beam;
lin  = cov;                                           % mascara de celdas validas

ownLin = nan(M, Nt);
idx3   = sub2ind([M nB Nt], mIdx(lin), bOwn(lin), kIdx(lin));
ownLin(lin) = Glin(idx3);

% Suma co-canal segun la sub-banda (regla unica)
sumSel = nan(M, Nt);
sA     = reshape(sumAll, M, Nt);
isInt  = lin & (sb == 0);
sumSel(isInt) = sA(isInt);
for c = 1:Delta
    sel = lin & (sb == c);
    if ~any(sel(:)), continue; end
    sC  = reshape(sumCol(:,c,:), M, Nt);
    sumSel(sel) = sC(sel);
end

% Interferencia intra = suma co-canal MENOS el propio haz (que es la portadora)
Iintra_lin_rel = max(sumSel - ownLin, 0);             % >= 0 (robustez numerica)

% Numero de haces co-canal (excluido el propio)
cntCol = accumarray(color, 1, [Delta 1]);
nCo    = nan(M, Nt);
nCo(isInt) = nB - 1;
for c = 1:Delta
    sel = lin & (sb == c);
    nCo(sel) = cntCol(c) - 1;
end

%% 3. Potencias absolutas
BuMHz_dB = 10*log10(Bu / 1e6);                        % [M x Nt]
Cpsd     = ctx.Cpsd_dBWMHz;

C_dBW = Cpsd + BuMHz_dB + 10*log10(ownLin);
N_dBW = -228.6 + 10*log10(cfg.radio.Tsys_K) + 10*log10(Bu);

I_intra_dBW = Cpsd + BuMHz_dB + 10*log10(Iintra_lin_rel);   % -Inf si no hay co-canal

% INTER-satelite: fraccion de coincidencia de sub-banda f
f = nan(M, Nt);
f(isInt) = 1;
f(lin & sb > 0) = 1 / Delta;
I_inter_dBW = ctx.Ipsd_inter_dBWMHz + BuMHz_dB + 10*log10(f);

%% 4. SINR
Nlin  = 10.^(N_dBW/10);
Ii    = 10.^(I_intra_dBW/10);   Ii(isnan(Ii)) = 0;
Ie    = 10.^(I_inter_dBW/10);   Ie(isnan(Ie)) = 0;
Clin  = 10.^(C_dBW/10);

SINR_lin = Clin ./ (Nlin + Ii + Ie);

FFR.SINR_dB      = 10*log10(SINR_lin);
FFR.CN_dB        = C_dBW - N_dBW;
FFR.CIR_intra_dB = 10*log10(ownLin ./ max(Iintra_lin_rel, realmin));
FFR.C_dBW        = C_dBW;
FFR.N_dBW        = N_dBW;
FFR.I_intra_dBW  = I_intra_dBW;
FFR.I_inter_dBW  = I_inter_dBW;
FFR.nCo          = nCo;
FFR.penalty_dB   = FFR.CN_dB - FFR.SINR_dB;

% Mascara de cobertura
flds = {'SINR_dB','CN_dB','CIR_intra_dB','C_dBW','N_dBW','I_intra_dBW', ...
        'I_inter_dBW','nCo','penalty_dB'};
for i = 1:numel(flds)
    FFR.(flds{i})(~cov) = NaN;
end

%% 5. Comprobaciones de cordura
% (a) La SINR nunca puede superar el C/N (la interferencia solo resta).
viol = FFR.SINR_dB - FFR.CN_dB;
if any(viol(cov) > 1e-6)
    warning('compute_sinr_ffr:sinr_gt_cn', ...
        'SINR > C/N en %d muestras (max +%.2e dB): revisar.', ...
        sum(viol(cov) > 1e-6), max(viol(cov)));
end
% (b) C/N debe ser INVARIANTE al ancho asignado (convencion de PSD constante):
%     C/N = Cpsd - 60 + 228.6 - 10*log10(Tsys) + Grel_propio, sin B_user.
%     Si esta identidad falla, se ha roto la convencion de potencia.
if any(cov(:))
    CN_id = ctx.Cpsd_dBWMHz(cov) - 60 + 228.6 - 10*log10(cfg.radio.Tsys_K) ...
            + 10*log10(ownLin(cov));
    dd = max(abs(FFR.CN_dB(cov) - CN_id));
    if dd > 1e-6
        warning('compute_sinr_ffr:cn_bandwidth', ...
            'C/N depende del ancho asignado (desvio %.2e dB): revisar la convencion de PSD.', dd);
    end
end
end
