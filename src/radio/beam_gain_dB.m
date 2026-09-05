function Grel_dB = beam_gain_dB(theta_off_deg, beamwidth3dB_deg, floor_dB)
%BEAM_GAIN_DB  Ganancia relativa del patron de haz Bessel (apertura circular).
%   Grel_dB = beam_gain_dB(theta_off_deg, beamwidth3dB_deg)
%   Grel_dB = beam_gain_dB(theta_off_deg, beamwidth3dB_deg, floor_dB)
%
%   Patron de apertura circular de 3GPP TR 38.811 sec. 6.4.1, parametrizado por
%   el ANCHO DE HAZ a -3 dB (no por el diametro de apertura):
%
%       u       = 1.61634 * sin(theta_off) / sin(theta_3dB/2)
%       Grel_dB = 20*log10(|2*J1(u)/u|)        (<= 0 dB, relativo a Gmax)
%
%   La constante 1.61634 es el PUNTO DE MEDIA POTENCIA del factor de apertura, es
%   decir la raiz de  |2*J1(u)/u| = 1/sqrt(2) = 0.707107  (verificado: en u = 1.61634
%   el factor vale 0.707107). Al normalizar u por sin(theta_3dB/2), en
%   theta_off = theta_3dB/2 sale u = 1.61634 y por tanto -3.0103 dB EXACTOS, sea cual
%   sea el ancho de haz (comprobado para HPBW = 0.5 / 1.0 / 1.5 / 1.7647 / 3.0 deg).
%   Ese es el sentido de la parametrizacion por ancho de haz a -3 dB.
%
%   CORRECCION (auditoria, hallazgo B9): la version anterior de este comentario decia
%   que 1.61634 era "el primer cero de d/du[2*J1(u)/u]". Es FALSO y conviene dejarlo
%   escrito porque es el tipo de dato que se pregunta en una defensa: los extremos
%   locales de |2*J1(u)/u| estan en u = 3.8317 (primer CERO del patron) y u = 5.1356
%   (pico del primer lobulo lateral, -17.57 dB). La formula y el valor numerico
%   siempre fueron correctos; lo unico erroneo era la justificacion.
%
%   Limite u->0 (boresight): 2*J1(u)/u -> 1, Grel_dB = 0.
%
%   floor_dB (OPCIONAL, por defecto -Inf): suelo de lobulos laterales relativo
%   al pico.  Grel_dB = max(Bessel, floor_dB). El patron Bessel ideal tiene nulos
%   profundos (-Inf) que infraestiman la interferencia fuera de eje; un suelo
%   realista (p.ej. -30 dB) la acota. El upgrade riguroso seria ITU-R S.1528
%   (trabajo futuro). Sin el 3er argumento, el comportamiento es identico al
%   anterior (callers de radioenlace no cambian).
%
%   La ganancia absoluta es  G = Gmax_dBi + Grel_dB , con Gmax INDEPENDIENTE.
%   Funcion reutilizable: haz del satelite (compute_link_budget, compute_interference).
%
%   Acepta entradas escalares o matriciales (theta_off_deg vectorizado).

if nargin < 3, floor_dB = -Inf; end

u = 1.61634 .* sind(theta_off_deg) ./ sind(beamwidth3dB_deg/2);

Grel_dB = zeros(size(u));                 % limite u->0 => 0 dB (boresight)
nz = abs(u) >= 1e-6;
Grel_dB(nz) = 20 * log10(abs(2 * besselj(1, u(nz)) ./ u(nz)));

Grel_dB = max(Grel_dB, floor_dB);         % suelo de lobulos laterales
Grel_dB(isnan(u)) = NaN;                  % NaN de entrada (sin cobertura) -> NaN
                                          % (sin esta linea, u=NaN caeria en la
                                          %  rama u->0 y devolveria 0 dB)
end
