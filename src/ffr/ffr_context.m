function ctx = ffr_context(cfg, sats, users, BL, R_ecef, G, S, L, INTF, verbose)
%FFR_CONTEXT  Contexto multihaz Earth-fixed: ganancias de haz y PSDs de referencia.
%   ctx = ffr_context(cfg, sats, users, BL, R_ecef, G, S, L, INTF)
%   ctx = ffr_context(..., INTF, verbose)   verbose = false: sin log (barridos
%         Monte Carlo / parfor, donde el log de cada punto solo es ruido). El
%         defecto es true, luego los llamadores existentes NO cambian.
%
%   PIEZA DE MOTOR compartida por TODOS los esquemas (reuso-1, reuso-Delta, FFR) y
%   por el barrido de alpha: precalcula UNA SOLA VEZ la fisica que NO depende del
%   esquema, para que ffr_allocate / compute_sinr_ffr / compute_kpis se limiten a
%   la logica de asignacion y los runners no contengan calculo.
%
%   MODELO DE HACES: **Earth-fixed** (decision del tutor). Las celdas son fijas en
%   el suelo (BL.cell_ecef, teselado hex de build_beam_layout centrado en
%   cfg.ground.point) y el satelite servidor apunta un haz a cada celda. Para el
%   usuario m en el instante k, servido por el satelite A:
%
%       psi(m,b,k) = angulo EN EL SATELITE entre  (A -> usuario m)  y
%                                                (A -> centro de celda b)
%       Grel(m,b,k) = beam_gain_dB(psi, HPBW, sidelobe_floor)      [dB, <= 0]
%
%   Es el angulo off-boresight EXACTO (no la aproximacion de plano tangente que
%   usa BL.cell_km, que solo define el teselado). Consecuencia geometrica que
%   importa a la defensa: al bajar la elevacion el cluster Earth-fixed se COMPRIME
%   angularmente visto desde el satelite -> los haces se solapan mas -> la
%   interferencia intra-satelite EMPEORA. Esta dependencia con la elevacion es
%   justo la que motiva la FFR ADAPTATIVA (H2, gancho en ffr_policy).
%   SIMPLIFICACION anotada: HPBW fijo (sin ensanchamiento del haz al apuntar fuera
%   de boresight ni celdas elipticas); igual que en compute_link_budget/E0.
%
%   Se REUTILIZAN sin recalcular:
%     L.FSPL_dB, L.Latm_clear_dB  (radioenlace del servidor)
%     INTF.I_dBW                  (interferencia INTER-satelite agregada, reuso-1)
%
%   CONVENCION DE POTENCIA: **PSD constante**. El satelite radia una densidad EIRP
%   fija (cfg.radio.EIRPdensity_dBWMHz) sobre el ancho que tenga asignado, luego la
%   portadora recibida en un ancho B_u es  Cpsd + 10*log10(B_u[MHz]) + Grel.  Asi
%   C/N NO depende del ancho asignado y la comparacion entre esquemas aisla el
%   compromiso real de la FFR (ancho de banda vs interferencia).
%   Alternativa (potencia total por haz constante -> boost de PSD al estrechar la
%   sub-banda, tipo SFR) = extension futura.        % <-- DECISION
%
%   Salidas (struct ctx):
%     .BL, .M, .Nt, .nBeams
%     .Grel        [M x nBeams x Nt] ganancia relativa de cada haz al usuario [dB]
%     .Glin        [M x nBeams x Nt] idem en lineal (evita repetir 10.^(x/10))
%     .beam        [M x Nt]  haz servidor (el de MAYOR ganancia hacia el usuario)
%     .theta_deg   [M x Nt]  offset angular del usuario a su centro de haz [deg]
%     .Cpsd_dBWMHz [M x Nt]  PSD de portadora de referencia (Grel = 0) [dBW/MHz]
%     .Ipsd_inter_dBWMHz [M x Nt] PSD interferente INTER-satelite (reuso-1) [dBW/MHz]
%     .elev_deg, .range_km [M x Nt]  geometria del servidor
%     .cov         [M x Nt]  logico: hay cobertura
%     .tvec                  vector de tiempos (informativo)
%     .SINR_ref_dB [M x Nt]  SINR de referencia reuso-1 (lo rellena ffr_allocate
%                            o el runner; se usa en el clasificador 'sinr')

if nargin < 10 || isempty(verbose), verbose = true; end

if ~isfield(BL,'cell_ecef')
    error('ffr_context:noCells', ...
        ['BL no trae la proyeccion Earth-fixed de celdas (BL.cell_ecef). ' ...
         'Necesita cfg.ground.point en build_beam_layout.']);
end

[M, Nt] = size(S.idx);
nB      = BL.nBeams;
HPBW    = cfg.radio.beamwidth3dB_deg;
fl      = cfg.radio.sidelobe_floor_dB;

ctx.BL     = BL;
ctx.M      = M;
ctx.Nt     = Nt;
ctx.nBeams = nB;
ctx.cov    = ~isnan(S.idx);

%% 1. Angulo off-boresight exacto usuario/haz visto desde el satelite servidor
psi = nan(M, nB, Nt);
C   = BL.cell_ecef;                         % [nB x 3] centros de celda (ECEF, km)

for k = 1:Nt
    Rk = R_ecef(:,:,k);                     % [N x 3] posiciones de satelites
    for m = 1:M
        A = S.idx(m,k);
        if isnan(A), continue; end
        rs = Rk(A,:);                       % satelite servidor
        du = users.ecef(m,:) - rs;          % satelite -> usuario
        db = C - rs;                        % satelite -> cada centro de celda
        nu = norm(du);
        nb = sqrt(sum(db.^2, 2));
        ca = (db * du.') ./ (nb * nu);      % [nB x 1] cosenos
        psi(m,:,k) = acosd( max(-1, min(1, ca)) ).';
    end
end

%% 2. Patron de haz (Bessel 38.811 + suelo de lobulos), vectorizado en el tensor
ctx.Grel = beam_gain_dB(psi, HPBW, fl);     % [M x nB x Nt], NaN donde no hay cobertura
ctx.Glin = 10.^(ctx.Grel/10);

%% 3. Haz SERVIDOR del usuario = el de mayor ganancia (menor offset angular)
%    Regla operativa: el terminal se asocia al haz que mejor lo ilumina.
[~, bidx]     = max(ctx.Grel, [], 2, 'omitnan');   % [M x 1 x Nt]
ctx.beam      = reshape(bidx, M, Nt);
[thmin, ~]    = min(psi, [], 2);
ctx.theta_deg = reshape(thmin, M, Nt);
ctx.beam(~ctx.cov)      = NaN;
ctx.theta_deg(~ctx.cov) = NaN;

%% 4. PSD de portadora de referencia (Grel = 0) y PSD interferente inter-satelite
%    Cpsd = EIRPdens + Grx_max - FSPL - Latm     [dBW/MHz]
%    (EIRPdens ya incluye la ganancia pico del haz; Grel se suma aparte).
%    El terminal apunta al satelite servidor -> Grx_max para la portadora Y para
%    la interferencia INTRA (mismo satelite): el VSAT no discrimina intra-satelite.
ctx.Cpsd_dBWMHz = cfg.radio.EIRPdensity_dBWMHz + cfg.radio.Grx_max_dBi ...
                  - L.FSPL_dB - L.Latm_clear_dB;

% INTER-satelite: INTF.I_dBW es potencia agregada sobre el ancho TOTAL en reuso-1
% -> se pasa a densidad para poder reescalarla al ancho realmente asignado.
ctx.Ipsd_inter_dBWMHz = INTF.I_dBW - 10*log10(cfg.radio.B_MHz);

%% 5. Geometria del servidor y tiempos
ctx.elev_deg = S.el;
ctx.range_km = S.range;
ctx.nVis     = S.nVis;
ctx.tvec     = [];
ctx.SINR_ref_dB = [];                       % lo rellena ffr_allocate / el runner

%% 6. Comprobaciones de cordura
% (a) El haz servidor debe ser casi siempre el de la celda geograficamente mas
%     cercana; se informa la fraccion de usuarios que caen en el haz central.
if any(ctx.cov(:))
    thv = ctx.theta_deg(ctx.cov);

    % (b) La rejilla no puede sobresalir del cluster de haces. El criterio (umbral
    %     = circunradio de celda s/sqrt(3), y su justificacion geometrica) vive en
    %     grid_outside_cluster, FUENTE UNICA compartida con run_e3b_minelev y
    %     run_convergence_nrings: el motor y los diagnosticos de los experimentos
    %     no pueden divergir. Es SOLO un aviso: no entra en ningun calculo ni
    %     cambia ningun KPI.
    %     `pctCirc` (offset maximo en % del circunradio) se reporta SIEMPRE, no solo
    %     al fallar: es la magnitud que separa limpiamente los puntos sanos (92-98%,
    %     saturan el vertice por debajo) de los sesgados (118-904%).
    %     El circunradio depende de la RETICULA (hex: s/sqrt(3); rect: s*sqrt(2)/2),
    %     asi que se le pasa BL.lattice. Los layouts historicos no traen ese campo y
    %     grid_outside_cluster asume 'hex' -> bit-identico.
    latt = 'hex';
    if isfield(BL,'lattice') && ~isempty(BL.lattice), latt = BL.lattice; end
    [isOut, pctCirc, thLim] = grid_outside_cluster(max(thv), BL.spacing_deg, latt);
    if strcmp(latt,'rect'), circLbl = 's*sqrt(2)/2'; else, circLbl = 's/sqrt(3)'; end

    if verbose
        fprintf('[ffr_context] M=%d usuarios | Nt=%d instantes | %d haces\n', M, Nt, nB);
        fprintf(['[ffr_context] offset angular al haz servidor: media %.3f deg, max %.3f deg ' ...
                 '(s=%.3f deg | circunradio %s=%.3f deg -> %.1f%%)\n'], ...
            mean(thv), max(thv), BL.spacing_deg, circLbl, thLim, pctCirc);
        fprintf('[ffr_context] elevacion del servidor: media %.1f deg, min %.1f deg, max %.1f deg\n', ...
            mean(ctx.elev_deg(ctx.cov)), min(ctx.elev_deg(ctx.cov)), max(ctx.elev_deg(ctx.cov)));
    end

    if isOut
        warning('ffr_context:offset', ...
            ['Hay usuarios FUERA de su celda (max offset %.3f deg = %.1f%% del ' ...
             'circunradio de celda %s = %.3f deg, con s = %.3f deg): la rejilla de ' ...
             'usuarios sobresale del cluster de haces.'], ...
            max(thv), pctCirc, latt, thLim, BL.spacing_deg);
    end
end
end
