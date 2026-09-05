function R_ecef = eci2ecef(cfg, R_eci, tvec)
%ECI2ECEF  Rota posiciones ECI -> ECEF por el angulo de rotacion terrestre.
%   R_eci, R_ecef : [N x 3 x Nt]
%
%   r_ecef = Rz(-theta) * r_eci ,  con  theta = theta_g0 + we*t  (GMST).
%   Trabajar en ECEF permite que el punto de tierra sea FIJO (rota con la
%   Tierra) mientras los satelites se mueven: es la referencia "desde el suelo".

we  = cfg.const.we;
th0 = deg2rad(cfg.prop.theta_g0);
[N,~,Nt] = size(R_eci);
R_ecef = zeros(N,3,Nt);

for k = 1:Nt
    th = th0 + we*tvec(k);
    c = cos(th);  s = sin(th);
    X = R_eci(:,1,k);  Y = R_eci(:,2,k);  Z = R_eci(:,3,k);
    R_ecef(:,1,k) =  c.*X + s.*Y;     % Rz(-theta)
    R_ecef(:,2,k) = -s.*X + c.*Y;
    R_ecef(:,3,k) =  Z;
end
end
