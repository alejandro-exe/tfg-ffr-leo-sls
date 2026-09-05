function BL = build_beam_layout(cfg)
%BUILD_BEAM_LAYOUT  Reticula de haces (hex) en el marco de antena del satelite.
%   BL = build_beam_layout(cfg)
%
%   Construye un layout hexagonal de haces (centro + anillos) con la SEPARACION
%   INTER-HAZ normativa del 3GPP TR 38.821 Sec. 6.1.1:
%
%       s_deg = rad2deg( sqrt(3) * sin( deg2rad(HPBW)/2 ) )
%
%   donde HPBW = cfg.radio.beamwidth3dB_deg. Esta es la separacion entre centros
%   de haces adyacentes en el dominio angular (aprox. de pequenos angulos de la
%   forma exacta asin(sqrt(3)*sin(HPBW/2))). El haz central esta en el nadir del
%   satelite (offset 0,0); los demas forman anillos hexagonales alrededor.
%
%   AVISO -- EL USO **INVERSO** DE ESTA RELACION NO ES INOCUO
%   ---------------------------------------------------------
%   Hacia adelante la aproximacion es excelente en el regimen para el que se
%   dedujo (angulo pequeno): con HPBW = 2.08 deg (perfil Ku) las dos formas dan
%   s = 1.801234 vs 1.801531 deg, un 0.0165 %. Pero el error CRECE con el ancho
%   de haz, y quien INVIERTA la relacion para deducir la HPBW a partir de una
%   separacion (o de una huella) dada lo arrastra amplificado:
%
%       inversa aproximada   HPBW = 2*asin( s[rad] / sqrt(3) )
%       inversa exacta       HPBW = 2*asin( sin(s) / sqrt(3) )
%
%   Ejemplo REAL del proyecto (config_oneweb: celda de 285 km desde 1200 km ->
%   s = 13.360218 deg): la inversa aproximada da 15.4740 deg y la exacta 15.3333
%   deg, un 0.91 %. Es el UNICO punto del trabajo donde la relacion se usa fuera
%   del regimen de angulo pequeno. NO invertirla a ciegas.
%
%   Y si alguna vez se pasa a la forma exacta, hay que cambiar LAS DOS A LA VEZ
%   (esta funcion y quien la invierta). Hoy se cancelan entre si, y por eso el
%   layout de OneWeb reproduce su huella declarada EXACTAMENTE
%   (BL.spacing_km = 285.000000 km). Cambiar solo la inversa DES-ANCLA la huella
%   (285 -> 282.33 km); cambiar solo esta funcion mueve la s de TODOS los
%   perfiles (Starlink 1.801234 -> 1.801531 deg) y rompe la reproducibilidad bit
%   a bit del corpus. La sensibilidad medida del cambio es <= 0.22 dB de C/I en el
%   peor usuario y ~0.002 dB sobre SINR_edge_p5: por eso se deja como esta.
%
%   NOTA (defensa): NO se usa la convencion "1.0 anchos de haz" (que daria
%   s_deg = HPBW). Con HPBW = 1.7647 deg -> s_deg ~ 1.528 deg (haces solapados,
%   como en 38.821). El apuntamiento entre haces es fijo (Earth-fixed cluster).
%
%   TAMANO del cluster (dos formas, excluyentes):
%     cfg.beams.nRings = R   -> anillos completos: nBeams = 3*R*(R+1) + 1
%                               (R=2 -> 19, R=3 -> 37, R=4 -> 61, R=5 -> 91)
%     cfg.beams.nRings = []  -> se usa cfg.beams.nBeams tal cual (posible anillo
%                               exterior incompleto). Es el caso de la CALIBRACION
%                               E0 (nBeams = 19, que ademas es anillo completo).
%
%   Salidas (struct BL):
%     BL.nBeams        - numero de haces
%     BL.spacing_deg   - separacion inter-haz s_deg [deg]
%     BL.offset_deg    - [nBeams x 2] offset angular (x,y) de cada haz respecto al
%                        haz central, en el plano angular del satelite [deg].
%                        Fila 1 = haz central (0,0).
%     BL.ring          - [nBeams x 1] anillo hex de cada haz (0 = central)
%     BL.ij            - [nBeams x 2] indices AXIALES (i,j) de la reticula hex
%                        (base e1 a 0 deg, e2 a 60 deg): p = i*e1 + j*e2.
%                        Los usa ffr_coloring para asignar sub-banda por haz.
%
%   PROYECCION EARTH-FIXED a tierra (celdas fijas en el suelo; los haces siguen a
%   su celda). Se calcula si cfg.ground tiene punto de referencia:
%     BL.h_ref_km      - altitud de referencia usada en la proyeccion [km]
%     BL.spacing_km    - separacion inter-celda proyectada [km]
%     BL.cell_km       - [nBeams x 2] centro de celda en coords locales (Este,Norte) [km]
%     BL.cell_lat/lon  - [nBeams x 1] centro de celda geodesico [deg]
%     BL.cell_ecef     - [nBeams x 3] centro de celda en ECEF [km]
%   La proyeccion es la del NADIR (plano tangente): d_km = h*tan(rho_deg), con
%   escalado RADIAL para que sea invariante a rotaciones. Define el teselado
%   Earth-fixed de celdas; los angulos reales usuario/haz vistos desde el
%   satelite los calcula ffr_context con la geometria exacta (no con esta aprox.).
%
%   Aportacion propia: generacion del cluster hexagonal y su proyeccion a celdas
%   Earth-fixed. La separacion sigue el estandar 3GPP TR 38.821 Sec. 6.1.1.

% --- Numero de haces: por anillos completos o por nBeams directo ---
nRings = [];
if isfield(cfg,'beams') && isfield(cfg.beams,'nRings'), nRings = cfg.beams.nRings; end
if ~isempty(nRings)
    nBeams = 3*nRings*(nRings+1) + 1;
else
    nBeams = cfg.beams.nBeams;
end
HPBW = cfg.radio.beamwidth3dB_deg;

% --- Topologia de la reticula: 'hex' (defecto) | 'rect' ---
% El modo 'rect' es un ANADIDO para el contraste de TOPOLOGIA de OneWeb (16 haces
% en 4x4). Si cfg.beams.lattice no existe se usa 'hex', luego TODOS los llamadores
% previos quedan BIT-IDENTICOS.
lattice = 'hex';
if isfield(cfg,'beams') && isfield(cfg.beams,'lattice') && ~isempty(cfg.beams.lattice)
    lattice = lower(cfg.beams.lattice);
end

% --- Separacion inter-haz normativa (3GPP TR 38.821 Sec. 6.1.1) ---
s_deg = rad2deg( sqrt(3) * sin( deg2rad(HPBW)/2 ) );

switch lattice
    case 'hex'
        % --- Generacion de centros hexagonales por anillos ---
        % Vectores base (separacion s_deg entre vecinos): e1 a 0 deg, e2 a 60 deg.
        % Los puntos son  p = a*e1 + b*e2  (a,b enteros).
        e1 = s_deg * [1, 0];
        e2 = s_deg * [cosd(60), sind(60)];

        % Enumeramos anillos hex crecientes hasta reunir nBeams haces.
        xy   = [0, 0];          % haz central
        ij   = [0, 0];          % indices axiales del haz central
        ring = 0;
        r    = 1;
        while size(xy,1) < nBeams
            % Los 6*r puntos del anillo r se recorren empezando en a=r,b=0 y girando.
            % Construccion estandar de anillo hexagonal en coordenadas axiales (a,b):
            aq = r; bq = 0;
            % 6 direcciones de avance sobre el borde del hexagono (axial):
            dirs = [ -1  1;  -1  0;   0 -1;   1 -1;   1  0;   0  1];
            for side = 1:6
                for step = 1:r
                    p = aq*e1 + bq*e2;
                    xy(end+1,:)   = p;          %#ok<AGROW>
                    ij(end+1,:)   = [aq, bq];   %#ok<AGROW>
                    ring(end+1,1) = r;          %#ok<AGROW>
                    if size(xy,1) >= nBeams, break; end
                    aq = aq + dirs(side,1);
                    bq = bq + dirs(side,2);
                end
                if size(xy,1) >= nBeams, break; end
            end
            r = r + 1;
        end

    case 'rect'
        % --- Reticula RECTANGULAR (contraste de topologia; OneWeb 4x4) ---
        % Base ortogonal de paso s_deg. Los indices (i,j) siguen siendo ENTEROS y
        % siguen cumpliendo  p = i*e1 + j*e2, que es todo lo que ffr_coloring
        % necesita: sus dos reglas (mod(i-j,3)+1 y mod(i,2)+2*mod(j,2)+1) son
        % coloreados VALIDOS tambien sobre reticula cuadrada (ningun par adyacente
        % comparte color). Lo que SI cambia es la distancia minima entre haces del
        % mismo color, que ffr_coloring mide y reporta:
        %     hex  Delta=3 -> sqrt(3)*s   |   rect Delta=3 -> sqrt(2)*s
        %     hex  Delta=4 -> 2*s         |   rect Delta=4 -> 2*s
        e1 = s_deg * [1, 0];
        e2 = s_deg * [0, 1];

        dims = [4 4];
        if isfield(cfg.beams,'rectDims') && ~isempty(cfg.beams.rectDims)
            dims = cfg.beams.rectDims;
        end
        if prod(dims) ~= nBeams
            error('build_beam_layout:rectDims', ...
                'cfg.beams.rectDims = [%d %d] da %d haces, pero nBeams = %d.', ...
                dims(1), dims(2), prod(dims), nBeams);
        end
        % Centros simetricos respecto al nadir: para un lado PAR no hay haz en (0,0),
        % asi que los indices van en semienteros escalados a enteros mediante un
        % desplazamiento comun. Se conserva la paridad de (i,j) -- que es lo que usa
        % el coloreado Delta=4 -- generando los indices como enteros 0..n-1 y
        % desplazando SOLO las coordenadas cartesianas.
        [ii, jj] = ndgrid(0:dims(1)-1, 0:dims(2)-1);
        ij = [ii(:), jj(:)];
        c  = (dims - 1)/2;                                   % centro geometrico
        xy = (ij(:,1)-c(1)).*e1 + (ij(:,2)-c(2)).*e2;
        % "Anillo" = anillo de Chebyshev respecto al centro (para diagnostico y para
        % que las funciones que agrupan por BL.ring sigan teniendo algo coherente).
        ring = max(abs(ij(:,1)-c(1)), abs(ij(:,2)-c(2)));
        ring = round(ring - min(ring));

        % El haz mas cercano al nadir se pone el PRIMERO, para conservar la
        % convencion "fila 1 = haz de referencia" del resto del motor.
        [~, i0] = min(sum(xy.^2, 2));
        ord = [i0; setdiff((1:size(xy,1)).', i0, 'stable')];
        xy = xy(ord,:);  ij = ij(ord,:);  ring = ring(ord);

    otherwise
        error('build_beam_layout:lattice', ...
            'cfg.beams.lattice = ''%s'' no soportado (usa ''hex'' o ''rect'').', lattice);
end

BL.nBeams      = nBeams;
BL.spacing_deg = s_deg;
BL.offset_deg  = xy(1:nBeams, :);
BL.ring        = ring(1:nBeams);
BL.ij          = ij(1:nBeams, :);
BL.lattice     = lattice;

% --- Comprobacion de cordura ---
if strcmp(lattice,'hex')
    % El haz central debe estar en (0,0); las distancias del primer anillo ~ s_deg.
    if any(BL.offset_deg(1,:) ~= 0)
        warning('build_beam_layout:center', 'El haz 1 no esta en el centro (0,0).');
    end
    d1 = sqrt(sum(BL.offset_deg(BL.ring==1,:).^2, 2));
    if ~isempty(d1) && max(abs(d1 - s_deg)) > 1e-9
        warning('build_beam_layout:ring1', ...
            'Anillo 1 a distancia %.4f (esperado s_deg=%.4f).', mean(d1), s_deg);
    end
    % Coherencia de los indices axiales con las posiciones cartesianas:
    xy_chk = BL.ij(:,1)*e1 + BL.ij(:,2)*e2;
    if max(abs(xy_chk(:) - BL.offset_deg(:))) > 1e-9
        warning('build_beam_layout:ij', 'BL.ij incoherente con BL.offset_deg.');
    end
else
    % En 'rect' el vecino mas proximo debe estar a exactamente s_deg.
    D = sqrt((BL.offset_deg(:,1)-BL.offset_deg(:,1).').^2 + ...
             (BL.offset_deg(:,2)-BL.offset_deg(:,2).').^2);
    D(1:nBeams+1:end) = Inf;
    if abs(min(D(:)) - s_deg) > 1e-9
        warning('build_beam_layout:rectSpacing', ...
            'Vecino mas proximo a %.4f deg (esperado s_deg=%.4f).', min(D(:)), s_deg);
    end
end

% --- Proyeccion Earth-fixed de las celdas al suelo (opcional) ---
if isfield(cfg,'ground') && isfield(cfg.ground,'point') && ~isempty(cfg.ground.point)
    h_ref = [];
    if isfield(cfg,'beams') && isfield(cfg.beams,'h_ref_km'), h_ref = cfg.beams.h_ref_km; end
    if isempty(h_ref), h_ref = cfg.constellations(1).h; end

    rho = sqrt(sum(BL.offset_deg.^2, 2));                 % offset radial [deg]
    sc  = ones(size(rho));                                % escala radial deg->km
    nz  = rho > 0;
    sc(nz) = h_ref * tand(rho(nz)) ./ rho(nz);            % km por grado (radial)
    cell_km = BL.offset_deg .* sc;                        % [nBeams x 2] (Este,Norte)

    lat0 = cfg.ground.point(1);  lon0 = cfg.ground.point(2);
    cell_lat = lat0 + cell_km(:,2) / 111.32;
    cell_lon = lon0 + cell_km(:,1) ./ (111.32 * cosd(cell_lat));

    BL.h_ref_km   = h_ref;
    BL.spacing_km = h_ref * tand(s_deg);
    BL.cell_km    = cell_km;
    BL.cell_lat   = cell_lat;
    BL.cell_lon   = cell_lon;
    BL.cell_ecef  = geodetic2ecef(cfg, cell_lat, cell_lon, cfg.ground.alt);
end
end
