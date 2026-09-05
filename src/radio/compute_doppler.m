function dop = compute_doppler(cfg, sats, S, R_ecef, V_eci, tvec)
%COMPUTE_DOPPLER  Corrimiento Doppler del satelite servidor sobre el usuario fijo.
%   Calcula el desplazamiento de frecuencia debido al movimiento orbital del
%   satelite servidor respecto a un punto FIJO en tierra (vista "desde el suelo").
%   Critico en LEO: las velocidades relativas son de varios km/s.
%
%   Aportacion propia: paso de velocidad ECI->ECEF (con termino de rotacion
%   terrestre), range rate y conversion a corrimiento Doppler.
%
%   Entradas:
%     cfg    -  configuracion (usa cfg.const.we, cfg.radio.freq_GHz, cfg.prop.theta_g0)
%     sats   -  constelacion (no usado directamente; se mantiene por simetria de API)
%     S      -  servidor: S.idx [M x Nt] indice del sat servidor (NaN sin cobertura)
%     R_ecef -  posiciones ECEF [N x 3 x Nt] [km]
%     V_eci  -  velocidades ECI [N x 3 x Nt] [km/s]  (2a salida de propagate)
%     tvec   -  vector de tiempos [1 x Nt] [s]
%
%   Salidas (struct dop, [M x Nt], NaN donde no hay cobertura):
%     df_Hz   - corrimiento Doppler [Hz]   (df = -rdot*f_c/c)
%     rdot_ms - range rate [m/s]
%
%   Convenio de signos:
%     rdot < 0  -> el satelite se ACERCA -> df > 0 (frecuencia recibida sube)
%     rdot = 0  -> acercamiento minimo (maxima elevacion) -> df = 0
%     rdot > 0  -> el satelite se ALEJA  -> df < 0

c    = 299792458;                 % m/s
f_c  = cfg.radio.freq_GHz * 1e9;  % Hz
we   = cfg.const.we;              % rad/s  rotacion terrestre
th0  = deg2rad(cfg.prop.theta_g0);

[M, Nt] = size(S.idx);
dop.df_Hz   = nan(M, Nt);
dop.rdot_ms = nan(M, Nt);

% Usuarios fijos en ECEF [km] (deterministas desde cfg; v_user = 0 en ECEF).
% El usuario es fijo en ECEF, asi que toda la velocidad relativa proviene del
% satelite (incluido el arrastre por la rotacion terrestre al pasar a ECEF).
users  = build_user_grid(cfg);
r_user = users.ecef;              % [M x 3] km

for m = 1:M
    ru = r_user(m, :);            % 1 x 3  (ECEF, km)
    for k = 1:Nt
        j = S.idx(m, k);
        if isnan(j), continue; end       % sin cobertura -> NaN

        th = th0 + we*tvec(k);
        cz = cos(th);  sz = sin(th);

        % --- Velocidad del satelite en ECEF (teorema de transporte) ---
        %   v_ecef = Rz(-theta)*v_eci  -  omega_e x r_ecef
        %   con omega_e = [0;0;we]  =>  omega_e x r = [-we*ry; we*rx; 0]
        v_eci = squeeze(V_eci(j, :, k)).';   % 3 x 1 km/s
        r_e   = squeeze(R_ecef(j, :, k)).';  % 3 x 1 km

        % Rz(-theta) * v_eci  (misma convencion que eci2ecef)
        vrot = [  cz*v_eci(1) + sz*v_eci(2);
                 -sz*v_eci(1) + cz*v_eci(2);
                  v_eci(3) ];
        omega_cross_r = [ -we*r_e(2);  we*r_e(1);  0 ];
        v_sat_ecef = vrot - omega_cross_r;   % 3 x 1 km/s

        % --- Range rate hacia el punto fijo (v_user = 0) ---
        r_rel = r_e - ru.';                  % 3 x 1 km  (sat - usuario)
        v_rel = v_sat_ecef;                  % 3 x 1 km/s
        rdot_kms = dot(r_rel, v_rel) / norm(r_rel);   % km/s

        rdot_ms = rdot_kms * 1e3;            % m/s
        dop.rdot_ms(m, k) = rdot_ms;
        dop.df_Hz(m, k)   = -rdot_ms * f_c / c;       % Hz
    end
end
end
