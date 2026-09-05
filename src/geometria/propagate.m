function [R_eci, V_eci] = propagate(cfg, sats, tvec)
%PROPAGATE  Posiciones (y opcionalmente velocidades) de los satelites en ECI.
%   R_eci : matriz [N x 3 x Nt]  posicion [km]   (N satelites, 3 coord, Nt instantes)
%   V_eci : matriz [N x 3 x Nt]  velocidad [km/s] (salida OPCIONAL, solo si se pide)
%
%   Orbita circular (e = 0): el argumento de latitud es  u = M0 + n*t .
%   Posicion en el plano orbital:  [a*cos(u) ; a*sin(u) ; 0] ,
%   y se lleva a ECI con   r_eci = Rz(RAAN) * Rx(inc) * r_perifocal .
%
%   Velocidad (derivada de la posicion, orbita circular):
%     v_perifocal = a*n*[-sin(u) ; cos(u) ; 0]  (km/s)
%   rotada con la MISMA transformacion Rz(RAAN)*Rx(inc) que la posicion.
%   Necesaria para el modulo Doppler (corrimiento por movimiento orbital, LEO).
%
%   SGP4: **NO IMPLEMENTADO**. Es TRABAJO FUTURO declarado en la memoria, no un
%   modo que se pueda activar por configuracion. Implementarlo exigiria (a) TLEs
%   reales de la constelacion -- que este TFG no usa, porque construye la Walker
%   sinteticamente desde cfg.constellations -- y (b) la funcion sgp4 de Aerospace
%   Toolbox, manteniendo el mismo formato de salida [N x 3 x Nt]. Ningun mecanismo
%   identificado en el TFG depende del propagador, asi
%   que no esta en la ruta critica. NO INTENTAR "activarlo" poniendo
%   cfg.prop.method = 'sgp4': no hay codigo detras, y por eso esta rutina aborta.

if ~strcmpi(cfg.prop.method,'kepler')
    error('propagate:notImplemented', ...
        ['cfg.prop.method = "%s" NO ESTA IMPLEMENTADO. El unico propagador del ' ...
         'simulador es "kepler" (orbita circular, sin J2 ni arrastre).\n' ...
         'SGP4 es TRABAJO FUTURO, no una opcion desactivada: no existe codigo que ' ...
         'lo respalde. Implementarlo requiere TLEs reales (el simulador construye ' ...
         'la constelacion Walker sinteticamente, no desde TLEs) y sgp4 de Aerospace ' ...
         'Toolbox, devolviendo [N x 3 x Nt]. Ver INFORME_TECNICO Sec.12.1 (L6) y ' ...
         'Sec.13.'], cfg.prop.method);
end

N = sats.N;  Nt = numel(tvec);
R_eci = zeros(N,3,Nt);
want_vel = nargout > 1;                   % solo calcula velocidad si se pide
if want_vel, V_eci = zeros(N,3,Nt); end

ci = cos(sats.inc);  si = sin(sats.inc);
cO = cos(sats.raan); sO = sin(sats.raan);

for k = 1:Nt
    u  = sats.M0 + sats.n * tvec(k);     % N x 1  argumento de latitud
    xp = sats.a .* cos(u);               % perifocal x
    yp = sats.a .* sin(u);               % perifocal y  (z = 0)

    % Rx(inc) * [xp; yp; 0] = [xp ; yp*ci ; yp*si]
    x1 = xp;       y1 = yp.*ci;   z1 = yp.*si;
    % Rz(RAAN) * [x1; y1; z1]
    X = cO.*x1 - sO.*y1;
    Y = sO.*x1 + cO.*y1;
    Z = z1;

    R_eci(:,:,k) = [X, Y, Z];

    if want_vel
        % v_perifocal = a*n*[-sin(u); cos(u); 0]  -> misma rotacion que la posicion
        vp_x = -sats.a .* sats.n .* sin(u);
        vp_y =  sats.a .* sats.n .* cos(u);
        % Rx(inc) * [vp_x; vp_y; 0]
        vx1 = vp_x;    vy1 = vp_y.*ci;   vz1 = vp_y.*si;
        % Rz(RAAN) * [vx1; vy1; vz1]
        VX = cO.*vx1 - sO.*vy1;
        VY = sO.*vx1 + cO.*vy1;
        VZ = vz1;
        V_eci(:,:,k) = [VX, VY, VZ];
    end
end
end
