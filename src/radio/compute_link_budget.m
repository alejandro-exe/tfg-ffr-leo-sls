function L = compute_link_budget(cfg, S)
%COMPUTE_LINK_BUDGET  Modulo de radioenlace downlink Ku (haz Earth-fixed).
%   Calcula el presupuesto de enlace para el satelite servidor en cada
%   instante temporal. La capa geometrica (S.range, S.el) es la entrada.
%
%   Aportacion propia: FSPL y patron de haz Bessel.
%   Toolbox (no aportacion): p618PropagationLosses para perdidas ITU-R P.618.
%
%   Entradas:
%     cfg  -  configuracion (requiere campos cfg.radio y cfg.geom.minElev)
%     S    -  satelite servidor: S.range [M x Nt] km, S.el [M x Nt] deg
%
%   Salidas (struct L, todos [M x Nt], NaN donde no hay cobertura):
%     FSPL_dB         - perdida en espacio libre [dB]
%     Gbf_dBi         - ganancia haz Bessel TX satelite [dBi]
%     Latm_clear_dB   - perdidas atmosfericas (50% excedencia, ~cielo claro) [dB]
%     Latm_rain_dB    - perdidas atmosfericas (p618_availability% excedencia) [dB]
%     CN_clear_dB     - C/N en cielo claro [dB]
%     CN_rain_dB      - C/N con lluvia [dB]
%     margin_clear_dB - margen de enlace cielo claro [dB]  (>0 = enlace viable)
%     margin_rain_dB  - margen de enlace con lluvia [dB]
%     Cshannon_Mbps   - capacidad Shannon (cielo claro, orientativo) [Mbps]
%     p618_used       - true si la toolbox P.618 estuvo disponible

%% 0. Coherencia de la configuracion (auditoria, hallazgo G3)
%  Punto de paso UNIVERSAL: todo pipeline del proyecto calcula balance de enlace
%  (E0, run_ffr_demo, E3, E3b, E4, E5, barridos), luego comprobar aqui el invariante
%  GT_dBK == Grx_max_dBi - 10log10(Tsys_K) cubre todos los caminos. Coste O(bloques
%  temporales), no O(elementos). Aborta si alguien sobrescribio uno de los dos sin
%  recalcular el otro (medido: GT_dBK = 20 dejaba Tsys obsoleto = 5 dB de error de
%  ruido, silencioso). Ver assert_cfg_coherent.
assert_cfg_coherent(cfg);

%% Constantes
[M, Nt] = size(S.el);

%% 1. FSPL (formula de Friis, programada a mano -> freespace_loss)
L.FSPL_dB = freespace_loss(S.range, cfg.radio.freq_GHz);

%% 2. Patron de haz Bessel (3GPP TR 38.811 sec. 6.4.1, apertura circular TX)
%    G(theta) = Gmax[dBi] + beam_gain_dB(theta_off, beamwidth3dB)  (funcion reutilizable)
%    Parametrizado por el ANCHO DE HAZ a -3 dB; Gmax INDEPENDIENTE de la apertura.
%    Con cellOffset_deg = 0 (usuario en el centro de celda) -> G = Gmax.
G_rel_dB  = beam_gain_dB(cfg.radio.cellOffset_deg, cfg.radio.beamwidth3dB_deg);
L.Gbf_dBi = (cfg.radio.Gmax_dBi + G_rel_dB) * ones(M, Nt);
L.Gbf_dBi(isnan(S.el)) = NaN;

%% 3. Perdidas atmosfericas ITU-R P.618-13 (via toolbox -> atm_loss_dB)
%    "Cielo claro" = excedencia 50%; "Lluvia" = cfg.radio.p618_availability%.
[L.Latm_clear_dB, L.Latm_rain_dB, L.p618_used] = atm_loss_dB(cfg, S.el);

%% 4. Balance C/N (formula derivada de Friis, programada a mano)
%    C/N[dB] = EIRPdens[dBW/MHz] + G/T[dB/K] + 168.6 - FSPL[dB] - Latm[dB]
%
%    Derivacion (unidades SI):
%      C/N = EIRP[dBW] + G/T - FSPL - Latm + 228.6 - B[dBHz]
%      EIRP[dBW] = EIRPdens + 10*log10(B_MHz)
%      B[dBHz]   = 10*log10(B_MHz*1e6) = 10*log10(B_MHz) + 60
%      -> B se cancela: C/N = EIRPdens + G/T - FSPL - Latm + 228.6 - 60
%                           = EIRPdens + G/T - FSPL - Latm + 168.6
%
%    La variacion del haz (delta_G = Gbf - Gmax) reduce el EIRP efectivo
%    para usuarios fuera del centro de celda (cellOffset_deg != 0).
%
%    *** NO "ARREGLES" ESTA CANCELACION *** (auditoria, hallazgo G2)
%    delta_G = (Gmax + Grel) - Gmax = Grel: cfg.radio.Gmax_dBi SE CANCELA EXACTO y
%    eso es CORRECTO, no un descuido. La densidad EIRP (cfg.radio.EIRPdensity_dBWMHz)
%    YA INCLUYE la ganancia de pico del haz -- es la definicion de EIRP --, luego lo
%    que entra en el enlace es SOLO la variacion RELATIVA del patron. Sumar Gmax otra
%    vez aqui, o anadirlo al Cpsd de ffr_context, seria DOBLE CONTABILIZACION.
%    Consecuencia concreta si alguien lo hiciera: en E3 la ganancia de pico se barre
%    como Gmax(theta) = Gmax_ref + 20log10(theta_ref/theta), asi que el C/N quedaria
%    inflado hasta +9.5 dB en el punto mas denso (theta = 0.5 deg) y la conclusion de
%    E3 -- que al densificar TODO lo que se mueve viene de la interferencia, porque el
%    C/N es invariante -- se vendria abajo sin que nada avisara.
delta_G      = L.Gbf_dBi - cfg.radio.Gmax_dBi;   % [dB], = Grel; Gmax se cancela (ver arriba)
base_offset  = cfg.radio.EIRPdensity_dBWMHz + cfg.radio.GT_dBK + 168.6;

L.CN_clear_dB = base_offset + delta_G - L.FSPL_dB - L.Latm_clear_dB;
L.CN_rain_dB  = base_offset + delta_G - L.FSPL_dB - L.Latm_rain_dB;

%% 5. Margen de enlace
L.margin_clear_dB = L.CN_clear_dB - cfg.radio.CNreq_dB;
L.margin_rain_dB  = L.CN_rain_dB  - cfg.radio.CNreq_dB;

%% 6. Capacidad Shannon orientativa (cielo claro)
%    C = B * log2(1 + SNR),  SNR_lineal = 10^(C/N_dB / 10)
L.Cshannon_Mbps = cfg.radio.B_MHz .* log2(1 + 10.^(L.CN_clear_dB / 10));

end
