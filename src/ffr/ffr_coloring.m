function [color, COL] = ffr_coloring(BL, Delta, verbose)
%FFR_COLORING  Coloreado de sub-bandas de la reticula hexagonal de haces.
%   color        = ffr_coloring(BL, Delta)
%   [color, COL] = ffr_coloring(BL, Delta, verbose)   (verbose=false: sin log)
%
%   PIEZA DE MOTOR de la FFR. Asigna a cada haz del layout (BL de
%   build_beam_layout) un COLOR 1..Delta, de modo que dos haces ADYACENTES nunca
%   compartan color. El color determina la sub-banda de borde que usa el haz, y
%   por tanto que haces son CO-CANAL entre si (los del mismo color). Es lo que
%   convierte la interferencia intra-satelite de reuso-1 en algo mitigable.
%
%   Se colorea sobre los indices AXIALES (i,j) de la reticula hex (BL.ij, base
%   e1 a 0 deg y e2 a 60 deg), no sobre coordenadas continuas: el coloreado es
%   exacto y sin ambiguedad numerica.
%
%     Delta = 1  ->  color = 1                     (reuso pleno; sin coloreado)
%     Delta = 3  ->  color = mod(i-j, 3) + 1       (patron clasico de reuso 3)
%     Delta = 4  ->  color = mod(i,2) + 2*mod(j,2) + 1
%
%   DISTANCIA MINIMA entre haces del MISMO color (en unidades de la separacion
%   inter-haz s), que es la que fija la reduccion de interferencia:
%     Delta = 3  ->  sqrt(3)*s ~ 1.732*s
%     Delta = 4  ->  2*s
%   La funcion la CALCULA y la IMPRIME (comprobacion de cordura del modulo).
%
%   Salidas:
%     color        - [nBeams x 1] color (sub-banda de borde) de cada haz, 1..Delta
%     COL (struct) - diagnostico:
%       .Delta          numero de colores
%       .count          [Delta x 1] haces por color
%       .dmin_s         distancia minima entre haces del mismo color, en unidades de s
%       .dmin_deg       idem en grados (marco angular del satelite)
%       .dmin_km        idem proyectada a tierra (si BL trae spacing_km)
%       .adjacentClash  true si algun par de haces ADYACENTES comparte color (fallo)
%
%   Aportacion propia del TFG (logica de reuso fraccionario).

if nargin < 3 || isempty(verbose), verbose = true; end

if ~isfield(BL,'ij')
    error('ffr_coloring:noIJ', ...
        'BL no trae indices axiales (BL.ij). Regenera el layout con build_beam_layout.');
end
i = BL.ij(:,1);  j = BL.ij(:,2);

switch Delta
    case 1
        color = ones(BL.nBeams, 1);
    case 3
        color = mod(i - j, 3) + 1;
    case 4
        color = mod(i, 2) + 2*mod(j, 2) + 1;
    otherwise
        error('ffr_coloring:Delta', ...
            'Delta = %d no soportado (usa 1, 3 o 4).', Delta);
end

%% --- Validacion: distancias entre haces del mismo color ---
s   = BL.spacing_deg;
D   = pdist2_local(BL.offset_deg);            % [nBeams x nBeams] distancias [deg]
sameColor = (color == color.');
D(1:BL.nBeams+1:end) = Inf;                   % ignorar la diagonal (mismo haz)

dmin_deg = min(D(sameColor));
if isempty(dmin_deg), dmin_deg = Inf; end     % Delta >= nBeams: un haz por color

% Pares ADYACENTES = a distancia s (tolerancia numerica).
% Con Delta = 1 (reuso pleno) TODOS los haces comparten color por definicion: no
% es un fallo del coloreado, es el peor caso espectral de referencia.
adj = abs(D - s) < 1e-6;
COL.adjacentClash = (Delta > 1) && any(adj(:) & sameColor(:));

COL.Delta    = Delta;
COL.count    = accumarray(color, 1, [Delta 1]);
COL.dmin_deg = dmin_deg;
COL.dmin_s   = dmin_deg / s;
COL.dmin_km  = NaN;
if isfield(BL,'spacing_km')
    COL.dmin_km = COL.dmin_s * BL.spacing_km;
end

if verbose
    fprintf('[ffr_coloring] Delta=%d | haces=%d | reparto por color: %s\n', ...
        Delta, BL.nBeams, mat2str(COL.count.'));
    if isnan(COL.dmin_km)
        fprintf('[ffr_coloring] distancia MINIMA entre haces del MISMO color: %.4f deg = %.3f*s\n', ...
            COL.dmin_deg, COL.dmin_s);
    else
        fprintf('[ffr_coloring] distancia MINIMA entre haces del MISMO color: %.4f deg = %.3f*s = %.2f km\n', ...
            COL.dmin_deg, COL.dmin_s, COL.dmin_km);
    end
end
if COL.adjacentClash
    warning('ffr_coloring:adjacent', ...
        'Hay haces ADYACENTES con el mismo color: el coloreado NO es valido.');
elseif verbose && Delta > 1
    fprintf('[ffr_coloring] OK: ningun par de haces adyacentes comparte color.\n');
end
end

% -------------------------------------------------------------------------
function D = pdist2_local(P)
%PDIST2_LOCAL  Matriz de distancias euclideas entre las filas de P (sin toolbox).
dx = P(:,1) - P(:,1).';
dy = P(:,2) - P(:,2).';
D  = sqrt(dx.^2 + dy.^2);
end
