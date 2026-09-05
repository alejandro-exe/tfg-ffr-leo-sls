function cfg = config_oneweb(layout)
%CONFIG_ONEWEB  Perfil de OneWeb como SISTEMA PROPIO (validacion cruzada de arquitecturas).
%   cfg = config_oneweb()          -> layout hexagonal de 19 haces (variante A, principal)
%   cfg = config_oneweb('hex')     -> idem
%   cfg = config_oneweb('rect')    -> reticula 4x4 de 16 haces (variante B, contraste)
%
%   PARTE DE config_default() y sobrescribe SOLO lo necesario, igual que config_calib
%   hace con el escenario de calibracion E0. config_default (Starlink Ku) queda INTACTO.
%
%   POR QUE ESTE PERFIL EXISTE
%   --------------------------
%   Starlink y OneWeb son arquitecturas radicalmente distintas: 550 vs 1200 km,
%   53 vs 87.9 deg de inclinacion, mascara de 25 vs 55 deg, y celdas de ~17 vs ~285 km.
%   Si las conclusiones del TFG (H1, mecanismo de compresion angular, criterio de
%   viabilidad) se reproducen en las DOS, dejan de ser propiedades del sistema
%   simulado y pasan a ser propiedades del REUSO MULTIHAZ EN LEO. Si alguna NO se
%   reproduce, el resultado es igual de valioso y acota el dominio de validez.
%
%   FUENTES (solicitudes ante la FCC y recomendaciones ITU-R; van en cada linea)
%     [Sched S]  = OneWeb Schedule S Technical Report
%     [TN p.NN]  = OneWeb Technical Narrative, pagina y seccion
%
%   DECISION DE MODELADO DEL LAYOUT (la clave del perfil)
%   ---------------------------------------------------------------------------------
%   OneWeb declara 16 haces Ku de usuario, ALTAMENTE ELIPTICOS y alineados norte-sur,
%   cubriendo una huella de 1140x1140 km [TN p.64 §A.1(iv), p.7 §A.2]. Eso plantea tres
%   problemas frente al motor, que asume celdas hexagonales y haces circulares:
%     (1) 16 haces no forman un cluster hexagonal completo (1+6+12 = 19);
%     (2) los haces son elipticos y beam_gain_dB modela apertura circular;
%     (3) la celda (~285 km) es un orden de magnitud mayor que la de Starlink.
%
%   Se ancla el layout en la HUELLA, no en Gmax:
%       s = 1140/4 = 285 km  ->  s_deg = atand(285/1200) = 13.3602 deg
%       HPBW = 2*asind(deg2rad(s_deg)/sqrt(3)) = 15.474 deg   (inversa de 38.821 §6.1.1)
%   Con 2 anillos hexagonales la separacion de celdas es 4*s = 4*285 = 1140 km, que es
%   la huella declarada, y las 19 celdas suman 1.34e6 km2 frente a los 1.30e6
%   declarados (+2.8%).
%     OJO al citar el tamano del cluster: 1140 km es 4*BL.spacing_km, NO la extension
%     PROYECTADA. La proyeccion deg->km de build_beam_layout es h*tan(rho), que no es
%     lineal, luego los centros del anillo 2 caen a 1200*tand(2*13.3602) = 604.07 km y
%     el cluster mide 1208.15 km de punta a punta (+6.0%). Las dos cifras son
%     correctas; lo que no vale es decir "1140 km de lado a lado".
%
%   AVISO SOBRE LA INVERSA (ver tambien la cabecera de build_beam_layout):
%   aqui la relacion normativa se usa en sentido INVERSO y fuera del regimen de angulo
%   pequeno para el que se dedujo. La inversa EXACTA, 2*asind(sind(s_deg)/sqrt(3)),
%   daria 15.3333 deg (0.91% menos). NO se corrige A PROPOSITO: la forma aproximada de
%   aqui y la de build_beam_layout SE CANCELAN, y por eso el layout reproduce la huella
%   declarada EXACTAMENTE (BL.spacing_deg = 13.3602184448 deg, BL.spacing_km =
%   285.000000 km). Corregir solo esta linea DES-ANCLA la huella (285 -> 282.33 km).
%   Sensibilidad medida del cambio: <= 0.22 dB de C/I en el peor usuario y ~0.002 dB
%   sobre SINR_edge_p5, frente a 1.90 dB hasta la frontera de veredicto mas ajustada
%   de E6 -> ningun veredicto puede cambiar.
%
%   SUPUESTO DECLARADO (discrepancia de 2 dB): la HPBW derivada de la huella implica
%   una directividad de ~22.4 dBi (10*log10(41253/HPBW^2) = 22.36, que es lo que
%   imprime run_oneweb_demo), frente a los 24.5 dBi declarados [Sched S p.23]. Se
%   ancla en la huella porque es la que fija la GEOMETRIA DE LA INTERFERENCIA, que es
%   lo que se estudia, y porque Gmax_dBi esta PROBADO que no entra en la fisica
%   (auditoria, hallazgo G2): la ganancia de pico ya esta contenida en
%   EIRPdensity_dBWMHz. Gmax_dBi se conserva solo como valor informativo.
%
%   POR QUE nRings = 2 ES EL MODELO CORRECTO, NO UNA APROXIMACION
%   ------------------------------------------------------------
%   Las celdas mas alla de la huella de un satelite OneWeb las sirven OTROS satelites
%   OneWeb, no mas haces del mismo. Y por coordinacion INTRA-operador esos vecinos NO
%   son co-canal con nosotros -- es exactamente el argumento con el que se justifica
%   cfg.interf.satPointing = 'nadir' (ver la cabecera de compute_interference). Por
%   tanto, modelar como interferencia intra SOLO los haces del satelite SERVIDOR es el
%   modelo fisicamente correcto. Subir nRings para "ganar guarda" INVENTARIA
%   interferencia co-canal que no existe, y ademas inflaria la SINR de forma
%   pesimista sin justificacion fisica.
%
%   Aportacion propia (perfil de sistema; la fisica es la del motor compartido).

if nargin < 1 || isempty(layout), layout = 'hex'; end

cfg = config_default();

%% ----- Constelacion: OneWeb -----------------------------------------------
% Walker STAR (polar), a diferencia del DELTA de Starlink.
cfg.constellations = struct( ...
    'name',    {'OneWeb'}, ...
    'pattern', {'star'},  ...   % polar                              [Sched S p.2]
    'h',       {1200},    ...   % km                                 [Sched S p.2]
    'inc',     {87.9},    ...   % deg                                [Sched S p.2]
    'T',       {720},     ...   % satelites                          [Sched S p.2]
    'P',       {18},      ...   % planos                             [Sched S p.2]
    'F',       {1},       ...
    'Om0',     {0},       ...
    'M0',      {0},       ...
    'minElev', {[]});           % [] a proposito: deja que mande cfg.geom.minElev y
                                % que el barrido de mascara (E3b equivalente) muerda.

%% ----- Geometria / visibilidad --------------------------------------------
% OneWeb no opera por debajo de 55 deg (45-50 permitido solo en bajas latitudes).
cfg.geom.minElev = 55;               % deg                           [TN p.13 §A.4]

%% ----- Radioenlace ---------------------------------------------------------
cfg.radio.freq_GHz = 12.0;           % GHz  misma banda Ku que Starlink (comparacion limpia)
cfg.radio.B_MHz    = 250;            % MHz  16 canales de 250 MHz     [TN p.7 §A.2]

% EIRP: -13.4 dBW/4 kHz -> +10*log10(1e6/4e3) = +23.98 dB -> 10.58 dBW/MHz
cfg.radio.EIRPdensity_dBWMHz = 10.58;   % dBW/MHz                     [TN p.20 §A.12]
cfg.radio.Gmax_dBi           = 24.5;    % dBi  INFORMATIVO, NO entra en la fisica (G2)  [Sched S p.23]

% Ancho de haz derivado de la HUELLA (ver cabecera). No se redondea: cualquier
% redondeo desplazaria s_km y con el la geometria de celdas.
s_km_target = 1140/4;                                   % 285 km      [TN p.64 §A.1(iv)]
s_deg       = atand(s_km_target / cfg.constellations(1).h);
cfg.radio.beamwidth3dB_deg = 2*asind( deg2rad(s_deg)/sqrt(3) );   % = 15.4740 deg

% ----- Terminal: EL MISMO que en Starlink (60 cm @ 12 GHz) -----------------
% OneWeb declara aperturas de 30-75 cm [TN]. Se usa 60 cm, identico al perfil de
% Starlink, para que la comparacion aisle la ARQUITECTURA ESPACIAL y no el terminal.
% A 12 GHz con 60 cm los valores son por tanto los mismos:
%   lambda = 0.024983 m -> D/lambda = 24.02 -> S.1428-1: G = 20log10(24.02)+7.7 = 35.31 dBi
%   Tsys   = 150 + 290*(10^(1.2/10)-1) = 242.3 K (TR 38.821 Tab. 6.1.1.1-3)
%   G/T    = 35.3 - 10log10(242.3) = 11.5 dB/K
cfg.radio.Grx_max_dBi = 35.3;        % dBi                            [ITU-R S.1428-1]
cfg.radio.GT_dBK      = 11.5;        % dB/K (MAESTRO en este perfil, como en config_default)
cfg.radio.Tsys_K      = 10^((cfg.radio.Grx_max_dBi - cfg.radio.GT_dBK)/10);   % DERIVADO

%% ----- Layout multihaz -----------------------------------------------------
switch lower(layout)
    case 'hex'
        % VARIANTE A (principal): 19 haces, 2 anillos hexagonales completos.
        cfg.beams.nRings = 2;
        cfg.beams.nBeams = [];
        cfg.beams.lattice = 'hex';
    case 'rect'
        % VARIANTE B (contraste de TOPOLOGIA): 4x4 = 16 haces exactos, retícula
        % cuadrada. ffr_coloring funciona tal cual sobre BL.ij cuadrado (Delta=3 y
        % Delta=4 siguen siendo coloreados validos), asi que el unico cambio de
        % motor es el modo 'rect' de build_beam_layout.
        cfg.beams.lattice = 'rect';
        cfg.beams.rectDims = [4 4];
        cfg.beams.nRings = [];
        cfg.beams.nBeams = 16;
    otherwise
        error('config_oneweb:layout', 'layout debe ser ''hex'' o ''rect'' (recibido ''%s'').', layout);
end
cfg.beams.h_ref_km = [];             % [] = altitud de la constelacion 1 (1200 km)

%% ----- Rejilla de usuarios -------------------------------------------------
% Con celdas de 285 km, el radio de 40 km de Starlink caeria DENTRO de una sola
% celda y no habria poblacion de borde. Se redimensiona para abarcar varias celdas
% conservando ~300 usuarios (el mismo tamano de muestra que Starlink):
%   radius 200 km -> guarda de 1.1 anillos con nRings=2 (maximo admisible 228 km)
%   step    20 km -> M ~ 314 usuarios;  step/s = 0.070
% NOTA sobre el criterio de convergencia (FASE A): el criterio se enuncia como
% FRACCION DEL TAMANO DE CELDA, step/s ~ 0.32-0.35. Aqui step/s = 0.070 es 5x MAS
% FINO que el convergido, luego sobra resolucion; el tamano de la rejilla lo fija
% el numero de usuarios, no la convergencia.
cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 200;
cfg.ground.step_km   = 20;

%% ----- FFR: anclas de la regla adaptativa acordes a la mascara -------------
% La regla alpha-elevacion interpola entre elMin_deg y elMax_deg. Con mascara de
% 55 deg, dejar el ancla inferior en 25 (valor de Starlink) haria que w(t) nunca
% bajase de 0.46 y la regla perderia casi todo su recorrido.
cfg.ffr.elMin_deg = cfg.geom.minElev;   % 55 deg
cfg.ffr.elMax_deg = 90;

%% ----- Comprobacion de coherencia -----------------------------------------
assert_cfg_coherent(cfg);

end
