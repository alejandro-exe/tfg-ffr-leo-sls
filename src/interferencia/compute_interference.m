function INTF = compute_interference(cfg, sats, users, G, S, L)
%COMPUTE_INTERFERENCE  Interferencia co-canal inter-satelite (reuso-1, peor caso).
%   INTF = compute_interference(cfg, sats, users, G, S, L)
%
%   Para cada usuario m e instante k con cobertura, el satelite servidor A
%   (S.idx) aporta la portadora C y TODOS los demas satelites visibles (reuso
%   pleno: comparten banda) aportan interferencia I. La SINR resultante es la
%   metrica de partida de los KPIs (percentil 5 de throughput, H1/H2).
%
%   APORTACION PROPIA del TFG (interferencia co-canal con doble discriminacion):
%     - patron del haz del satelite interferente (Bessel 38.811 + suelo de lobulos
%       cfg.radio.sidelobe_floor_dB), y
%     - patron del terminal RX en tierra (ITU-R S.1428, term_gain_dB; apunta al
%       servidor, los interferentes llegan fuera de eje y se atenuan).
%
%   SUPUESTOS de esta version (anotar en la defensa):
%     - Reuso-1 (cfg.interf.reuse = 1): peor caso espectral. La FFR (coloreado de
%       sub-bandas) reducira el conjunto de interferentes co-canal -> trabajo del
%       modulo de FFR.
%     - Interferencia INTER-satelite unicamente. La intra-satelite (multihaz del
%       mismo satelite) se trata en el modulo de FFR.
%     - Apuntamiento del interferente cfg.interf.satPointing:
%         'nadir'     (nominal)   -> haz centrado en su sub-satellite point.
%         'boresight' (cota sup.) -> haz apuntado al usuario (Gmax): peor caso.
%     - Patron Bessel del satelite con SUELO de lobulos (el ideal sin suelo tiene
%       nulos -Inf que infraestiman I). Upgrade riguroso: ITU-R S.1528 (futuro).
%     - Atmosfera de cielo claro (excedencia 50%) para C e I (coherente con L).
%
%   Entradas:
%     cfg, sats, users - configuracion, constelacion, usuarios
%     G   - geometria completa: G.el, G.az, G.range, G.vis  [M x N x Nt]
%     S   - servidor: S.idx, S.el  [M x Nt]
%     L   - presupuesto de enlace del servidor (usa L.FSPL_dB, L.Latm_clear_dB)
%
%   Salidas (struct INTF, [M x Nt], NaN donde no hay cobertura):
%     SINR_dB - relacion senal a ruido + interferencia [dB]
%     I_dBW   - potencia interferente agregada [dBW] (-Inf si nInterf = 0)
%     CN_dB   - C/N del servidor [dB] (= C - N; debe coincidir con L.CN_clear_dB)
%     nInterf - numero de interferentes co-canal [adim]

[M, Nt] = size(S.idx);
N       = sats.N;
Re      = cfg.const.Re;
EIRP_const = cfg.radio.EIRPdensity_dBWMHz + 10*log10(cfg.radio.B_MHz);  % EIRP por haz [dBW]
Grx_max    = cfg.radio.Grx_max_dBi;

% Ruido termico (coherente con compute_link_budget: Tsys derivado de G/T)
%   N[dBW] = -228.6 + 10*log10(Tsys) + 10*log10(B_Hz)
N_dBW = -228.6 + 10*log10(cfg.radio.Tsys_K) + 10*log10(radio_band_Hz(cfg));
N_lin = 10^(N_dBW/10);

% Precalculo vectorizado (una sola construccion de tabla P.618) para TODOS los
% satelites: FSPL y atmosfera de cielo claro vs su geometria.
FSPL_all          = freespace_loss(G.range, cfg.radio.freq_GHz);     % [M x N x Nt]
[Latm_clear_all]  = atm_loss_dB(cfg, G.el);                          % [M x N x Nt]

INTF.SINR_dB = nan(M, Nt);
INTF.I_dBW   = nan(M, Nt);
INTF.CN_dB   = nan(M, Nt);
INTF.nInterf = nan(M, Nt);

for m = 1:M
    for k = 1:Nt
        A = S.idx(m, k);
        if isnan(A), continue; end            % sin cobertura -> NaN

        % --- Portadora del servidor (reutiliza L: garantiza CN_dB == L.CN_clear) ---
        C_dBW = EIRP_const + Grx_max - L.FSPL_dB(m,k) - L.Latm_clear_dB(m,k);
        C_lin = 10^(C_dBW/10);
        INTF.CN_dB(m,k) = C_dBW - N_dBW;

        % --- Conjunto de interferentes co-canal (reuso-1): visibles, != A ---
        vis = squeeze(G.vis(m,:,k));          % 1 x N logico
        Bset = find(vis);  Bset(Bset==A) = [];
        INTF.nInterf(m,k) = numel(Bset);

        if isempty(Bset)
            INTF.I_dBW(m,k)   = -Inf;          % sin interferencia
            INTF.SINR_dB(m,k) = C_dBW - N_dBW; % SINR = C/N exactamente
            continue;
        end

        % reshape a fila [1 x nB] de forma robusta (squeeze+.' da columna y rompe
        % el broadcasting cuando hay >=2 interferentes).
        el_B   = reshape(G.el(m,Bset,k),         1, []);   % elevacion de cada interferente
        az_B   = reshape(G.az(m,Bset,k),         1, []);   % azimut
        h_B    = reshape(sats.h(Bset),           1, []);   % altitud
        FSPL_B = reshape(FSPL_all(m,Bset,k),     1, []);
        Latm_B = reshape(Latm_clear_all(m,Bset,k), 1, []);

        % (1) Patron del haz del interferente B segun su APUNTAMIENTO:
        %     'nadir'     (nominal): el haz apunta a su sub-satellite point; el
        %                 off-boresight hacia el usuario es el angulo nadir
        %                 eta = asin( Re/(Re+h) * cos(el) ).  Suelo de lobulos.
        %     'boresight' (cota sup.): el haz apunta AL usuario -> ganancia pico
        %                 (Grel = 0), maximiza la interferencia (peor caso real).
        switch lower(cfg.interf.satPointing)
            case 'boresight'
                Gsat_rel = zeros(1, numel(Bset));            % Gmax hacia el usuario
            otherwise   % 'nadir'
                eta = asind( (Re ./ (Re + h_B)) .* cosd(el_B) );
                Gsat_rel = beam_gain_dB(eta, cfg.radio.beamwidth3dB_deg, ...
                                        cfg.radio.sidelobe_floor_dB);
        end
        EIRP_B = EIRP_const + Gsat_rel;

        % (2) Patron del terminal RX (ITU-R S.1428): el terminal apunta al servidor
        %     A; cada interferente llega con separacion angular psi (vista desde el
        %     usuario).  term_gain_dB devuelve dBi ABSOLUTOS (psi=0 -> Grx_max).
        el_A = S.el(m,k);                      % = G.el(m,A,k)
        az_A = G.az(m,A,k);                    % azimut del servidor (S no lo trae)
        carg = sind(el_A).*sind(el_B) + cosd(el_A).*cosd(el_B).*cosd(az_A - az_B);
        carg = max(-1, min(1, carg));          % robustez numerica para acosd
        psi  = acosd(carg);
        Grx_B = term_gain_dB(cfg, psi);

        % (3) Potencia interferente de cada B y suma en lineal
        I_B_dBW = EIRP_B + Grx_B - FSPL_B - Latm_B;     % 1 x nB
        I_lin   = sum(10.^(I_B_dBW/10));

        INTF.I_dBW(m,k)   = 10*log10(I_lin);
        INTF.SINR_dB(m,k) = 10*log10(C_lin / (N_lin + I_lin));
    end
end

% --- Validacion interna: la SINR nunca puede superar el C/N ---
viol = INTF.SINR_dB - INTF.CN_dB;
tol  = 1e-6;
if any(viol(:) > tol)
    warning('compute_interference:sinr_gt_cn', ...
        'SINR > C/N en %d instantes (max +%.2e dB): revisar.', ...
        sum(viol(:) > tol), max(viol(:)));
end
end
