function CIR = compute_beam_cir(cfg, BL, linkC_dBW, N_dBW, I_inter_dBW)
%COMPUTE_BEAM_CIR  CIR/SINR intra-satelite del HAZ CENTRAL (metodologia 38.821).
%   CIR = compute_beam_cir(cfg, BL)
%   CIR = compute_beam_cir(cfg, BL, linkC_dBW, N_dBW, I_inter_dBW)
%
%   Nucleo de la CALIBRACION E0. Distribuye usuarios uniformemente en el haz
%   CENTRAL de un cluster multihaz (BL, de build_beam_layout) y calcula, por
%   usuario, la relacion Portadora/Interferencia co-canal INTRA-satelite:
%
%       C  = ganancia del haz central hacia el usuario (off-axis = |rho_u|)
%       I  = suma de ganancias de los OTROS haces del mismo satelite (reuso-1)
%       CIR = C / I
%
%   Por que la interferencia intra-satelite es la DOMINANTE (defensa): todos los
%   haces provienen del MISMO satelite, luego el terminal VSAT (que apunta al
%   satelite) los ve a todos EN EJE -> Grx_max comun a C e I y NO discrimina.
%   Ademas FSPL y atmosfera son comunes (mismo satelite) -> se cancelan y la CIR
%   queda como puro cociente de PATRONES de haz (coherente con 3GPP TR 38.821).
%
%   Distribucion de usuarios: uniforme en area dentro de un disco de radio
%   R_cell = cfg.beams.cellFrac * BL.spacing_deg (crossover con el haz vecino en
%   s/2). rho_u es el offset angular del usuario respecto al centro del haz.
%
%   SINR (opcional): si se pasan linkC_dBW (portadora del servidor en dBW en el
%   centro de celda), N_dBW (ruido termico) e I_inter_dBW (interferencia
%   inter-satelite agregada, o -Inf), se calcula tambien la SINR total del haz
%   central anadiendo ruido e interferencia inter-satelite:
%       C_u   = linkC_dBW + Grel_central(rho_u)
%       I_u   = I_intra(rho_u) [reescalado a dBW] + I_inter
%       SINR  = C_u / (N + I_u)
%
%   Salidas (struct CIR, vectores [nUsers x 1]):
%     CIR.rho_deg   - offset angular de cada usuario respecto al centro del haz
%     CIR.CIR_dB    - CIR intra-satelite (solo patrones) [dB]
%     CIR.SINR_dB   - SINR total (si se dieron los argumentos de enlace) [dB]
%     CIR.stats     - struct con media/mediana/p5/p50/p95 de CIR y SINR

nU    = cfg.beams.nUsers;
HPBW  = cfg.radio.beamwidth3dB_deg;
floorlobe = cfg.radio.sidelobe_floor_dB;
Rcell = cfg.beams.cellFrac * BL.spacing_deg;

% --- Distribucion uniforme en area dentro del disco de radio Rcell ---
% (semilla fija para reproducibilidad de la CDF de calibracion)
rng(0);
rr  = Rcell * sqrt(rand(nU,1));            % r ~ sqrt(U) -> uniforme en area
th  = 2*pi*rand(nU,1);
pu  = [rr.*cos(th), rr.*sin(th)];          % posicion angular [nU x 2] (deg)

% Centro del haz servidor = haz central del cluster (fila 1 de BL.offset_deg)
c_serv  = BL.offset_deg(1,:);              % (0,0)
c_intf  = BL.offset_deg(2:end,:);          % 18 haces interferentes

% --- Ganancias relativas de patron (Bessel + suelo) ---
% Portadora: off-axis del usuario respecto al centro del haz servidor
rho_serv = sqrt(sum((pu - c_serv).^2, 2));           % [nU x 1]
Grel_C   = beam_gain_dB(rho_serv, HPBW, floorlobe);  % [nU x 1] (<= 0 dB)

% Interferentes intra-sat: off-axis del usuario respecto a cada haz vecino
CIR_dB = nan(nU,1);
I_intra_rel_dB = nan(nU,1);                % Sum de ganancias relativas [dB]
for u = 1:nU
    doff = sqrt(sum((pu(u,:) - c_intf).^2, 2));       % [18 x 1] off-axis a cada haz
    Grel_I = beam_gain_dB(doff, HPBW, floorlobe);     % [18 x 1] (<= 0 dB)
    I_lin  = sum(10.^(Grel_I/10));
    I_intra_rel_dB(u) = 10*log10(I_lin);
    CIR_dB(u) = Grel_C(u) - I_intra_rel_dB(u);
end

CIR.rho_deg = rho_serv;
CIR.CIR_dB  = CIR_dB;

% --- SINR total opcional (anade ruido termico e interferencia inter-satelite) ---
CIR.SINR_dB = [];
if nargin >= 4 && ~isempty(linkC_dBW) && ~isempty(N_dBW)
    if nargin < 5 || isempty(I_inter_dBW), I_inter_dBW = -Inf; end
    C_u_dBW = linkC_dBW + Grel_C;                      % portadora por usuario
    % Interferencia intra en potencia absoluta: misma referencia que la portadora
    % pico (linkC_dBW ya es C en el centro de celda con Grel=0), luego
    %   I_intra_dBW(u) = linkC_dBW + I_intra_rel_dB(u)
    I_intra_dBW = linkC_dBW + I_intra_rel_dB;          % [nU x 1]
    I_tot_lin   = 10.^(I_intra_dBW/10) + 10.^(I_inter_dBW/10);
    N_lin       = 10.^(N_dBW/10);
    CIR.SINR_dB = 10*log10( 10.^(C_u_dBW/10) ./ (N_lin + I_tot_lin) );
end

% --- Estadisticos ---
CIR.stats.CIR_mean   = mean(CIR.CIR_dB);
CIR.stats.CIR_median = median(CIR.CIR_dB);
CIR.stats.CIR_p5     = prctile(CIR.CIR_dB, 5);
CIR.stats.CIR_p50    = prctile(CIR.CIR_dB, 50);
CIR.stats.CIR_p95    = prctile(CIR.CIR_dB, 95);
if ~isempty(CIR.SINR_dB)
    CIR.stats.SINR_mean   = mean(CIR.SINR_dB);
    CIR.stats.SINR_median = median(CIR.SINR_dB);
    CIR.stats.SINR_p5     = prctile(CIR.SINR_dB, 5);
    CIR.stats.SINR_p50    = prctile(CIR.SINR_dB, 50);
    CIR.stats.SINR_p95    = prctile(CIR.SINR_dB, 95);
end
end
