function cfg = config_default()
%CONFIG_DEFAULT  Configuracion central del simulador SLS LEO-FFR (capa geometrica).
%   Todos los parametros del experimento viven aqui. Cambiar de Starlink a
%   OneWeb, de intra- a inter-constelacion, o de punto fijo a rejilla, es solo
%   editar esta funcion. Las decisiones pendientes del tutor estan marcadas
%   con  % <-- DECISION
%
%   Unidades: distancias en km, angulos en grados, tiempo en segundos.

cfg = struct();

% ----- Constantes fisicas -----
cfg.const.mu = 398600.4418;        % km^3/s^2  GM de la Tierra
cfg.const.Re = 6378.137;           % km        radio ecuatorial WGS-84
cfg.const.e2 = 6.69437999014e-3;   % -         excentricidad^2 WGS-84
cfg.const.we = 7.2921159e-5;       % rad/s     velocidad de rotacion terrestre

% ----- Ventana temporal -----
cfg.time.t0       = 0;             % s
cfg.time.dt       = 10;            % s   paso temporal
cfg.time.duration = 2*3600;        % s   2 h (cubre varios pasos de satelite)

% ----- Constelaciones (lista; 1 = intra, >1 = inter-constelacion) -----  % <-- DECISION (saturacion)
% Notacion Walker  i:T/P/F .  pattern: 'delta' (Starlink) | 'star' (OneWeb, polar)
cfg.constellations = struct( ...
    'name',    {'Starlink-S1'}, ...
    'pattern', {'delta'}, ...
    'h',       {550},   ...   % km  altitud
    'inc',     {53},    ...   % deg inclinacion
    'T',       {1584},  ...   % total de satelites
    'P',       {72},    ...   % planos
    'F',       {1},     ...   % factor de fase de Walker (0..P-1)
    'Om0',     {0},     ...   % deg RAAN del primer plano
    'M0',      {0},     ...   % deg anomalia media del primer satelite
    'minElev', {[]});         % deg elevacion minima de OPERACION de ESTA constelacion
                              % [] = usar el global cfg.geom.minElev (ver build_constellation)
% IMPORTANTE: dejar minElev = [] aqui, NO un numero. Si se fijara, anularia el barrido
% de cfg.geom.minElev de E3b (la mascara por constelacion manda sobre la global). El
% campo existe aunque este vacio para que ampliar el array de constelaciones no falle
% con "Subscripted assignment between dissimilar structures".

% Para anadir OneWeb (caso inter-constelacion), usa un array de structs. OneWeb declara
% una elevacion minima de 55 deg, MUY distinta de los 25 de Starlink, y por eso la
% mascara es por constelacion (si no, entrarian como interferentes satelites ajenos a
% elevaciones a las que no transmiten):
% cfg.constellations(2) = struct('name','OneWeb','pattern','star', ...
%     'h',1200,'inc',87.9,'T',720,'P',18,'F',1,'Om0',0,'M0',0,'minElev',55);

% ----- Segmento terreno: rejilla de usuarios FIJA en ECEF -----  % <-- DECISION (alcance)
cfg.ground.mode  = 'point';        % 'point' | 'grid' | 'list' | 'localgrid'
cfg.ground.point = [40.33, -3.77]; % [lat lon] deg (Leganes / UC3M aprox.)
cfg.ground.alt   = 0.670;          % km  altura sobre el elipsoide
% Para mode = 'grid':
cfg.ground.latlim = [35 44];       % deg
cfg.ground.lonlim = [-10 4];       % deg
cfg.ground.dlat   = 1.0;           % deg resolucion en latitud
cfg.ground.dlon   = 1.0;           % deg resolucion en longitud
% Para mode = 'localgrid' (rejilla LOCAL alrededor de cfg.ground.point; la que
% necesita la FFR para tener usuarios de CENTRO y de BORDE de celda):
cfg.ground.radius_km = 40;         % km  radio de la rejilla local              % <-- DECISION
cfg.ground.step_km   = 3;          % km  resolucion de la rejilla local         % <-- DECISION
% AVISO de coste: la geometria es M*N*Nt. Con 'localgrid' (M ~ cientos) hay que
% reducir la constelacion y/o usar cfg.time.dt = 60 s o una ventana corta.

% ----- Geometria / visibilidad -----
% Elevacion minima GLOBAL (por defecto y suelo de la tabla P.618). Cada constelacion
% puede llevar la SUYA en cfg.constellations(c).minElev (ver build_constellation):
% OneWeb, por ejemplo, no opera por debajo de 55 deg. Si una constelacion no declara
% minElev, se le aplica este valor.
cfg.geom.minElev = 25;             % deg elevacion minima de servicio nominal   [SpaceX FCC]
cfg.geom.serving = 'maxElev';      % satelite servidor: 'maxElev' | 'minRange'
% Constelaciones ELEGIBLES como servidoras (indices de cfg.constellations).
% [] = todas (comportamiento historico, valido con UNA constelacion).
% En MULTI-OPERADOR hay que fijarlo (E4 usa 1): si no, el usuario se asociaria al
% mejor satelite de CUALQUIER operador -- se conectaria a la competencia -- y ademas
% serviria celdas proyectadas para la altitud de la constelacion 1. Ver
% associate_serving.
cfg.geom.serving_cid = [];                                                  % <-- DECISION

% ----- Propagacion -----
% UNICO METODO IMPLEMENTADO: 'kepler' (orbita circular, sin J2 ni arrastre).
% SGP4 NO ESTA IMPLEMENTADO y NO se ofrece como opcion: figura como TRABAJO FUTURO
% en la memoria (requiere TLEs reales + sgp4 de Aerospace Toolbox). propagate.m
% aborta con un error explicito ante cualquier otro valor. Se retiro 'sgp4' de
% esta lista a proposito: aparecia como si fuera seleccionable y no lo es.
cfg.prop.method   = 'kepler';      % 'kepler'  (unico valor admitido)
cfg.prop.theta_g0 = 0;             % deg  GMST en t0 (0 = simplificacion)

% ----- Radioenlace (downlink, banda Ku, haz Earth-fixed) -----
% PERFIL Ku REAL. Los valores ya NO son placeholders: proceden de las solicitudes
% presentadas por los operadores ante la FCC y de las recomendaciones ITU-R. La fuente
% va en el propio comentario (sustituye al antiguo "% <-- DECISION").
cfg.radio.freq_GHz           = 12.0;     % GHz  Ku DL, canales 10.7-12.7 GHz          [SpaceX FCC]
cfg.radio.B_MHz              = 250;      % MHz  ancho de canal de bajada (MAESTRO)    [SpaceX FCC]
% NO existe cfg.radio.B_Hz: era un campo DERIVADO que quedaba petrificado si alguien
% sobrescribia B_MHz despues (p.ej. via pt.cfg_over -> merge_cfg), y B_Hz alimenta el
% RUIDO TERMICO. El ancho en Hz se deriva SIEMPRE en el punto de uso con
% radio_band_Hz(cfg). Ver la cabecera de esa funcion (auditoria, hallazgo G3).
% Densidad EIRP REAL: maxima EIRP oblicua (slant) declarada para 540 km, -15.6 dBW/4 kHz.
% Conversion: -15.6 + 10*log10(1e6/4e3) = -15.6 + 23.98 = 8.38 dBW/MHz.
% (valores anteriores: 34 -> 12 dBW/MHz, ambos estimados; este es el de la solicitud)
cfg.radio.EIRPdensity_dBWMHz = 8.38;     % dBW/MHz  max slant @540 km (-15.6 dBW/4kHz)
                                         %          [SAT-MOD-20200417-00037 Tech Attach, Tab. A.3.1-1]
% Gmax_dBi es INFORMATIVO: NO ENTRA EN LA FISICA (auditoria, hallazgo G2).
% La ganancia de pico del haz ya esta CONTENIDA en EIRPdensity_dBWMHz (por definicion
% de EIRP), y sumarla aparte seria doble contabilizacion. De hecho se cancela exacto
% en compute_link_budget (delta_G = Gbf - Gmax = Grel) y NO aparece en ffr_context.
% Se conserva porque documenta la directividad del haz y lo imprimen los runners.
% CONSECUENCIA PRACTICA: al barrer el ancho de haz (E3/E3b), el acoplamiento de
% directividad 20log10(theta_ref/theta) hay que aplicarlo a EIRPdensity_dBWMHz, NO a
% Gmax_dBi -- un override de Gmax_dBi es un NO-OP silencioso.
cfg.radio.Gmax_dBi           = 38;       % dBi  ganancia pico haz satelite TX (INFORMATIVO, no entra en la fisica)  % <-- DECISION
% Ancho de haz REAL: dimensionado para celda de ~20 km de diametro desde 550 km.
% Con la separacion normativa s = rad2deg(sqrt(3)*sin(HPBW/2)) de 38.821 Sec. 6.1.1,
% HPBW = 2.08 deg -> s = 1.8013 deg -> s_km = 550*tand(s) = 17.30 km entre centros
% -> celda hexagonal de diametro equivalente en area ~18.2 km. (antes: 1.5 deg estimado)
cfg.radio.beamwidth3dB_deg   = 2.08;     % deg  ancho de haz a -3 dB -> celda ~20 km desde 550 km
cfg.radio.aperture_r_m       = 1.0;      % m    radio de apertura -- NO USADO (el patron se parametriza por beamwidth3dB_deg)
% G/T REAL del terminal: G/T = Grx_max - 10*log10(Tsys), con Grx_max = 35.3 dBi
% (terminal de 60 cm, abajo) y Tsys = 242.3 K de la cadena de ruido de
% 3GPP TR 38.821 Tabla 6.1.1.1-3 (Tant = 150 K, NF = 1.2 dB):
%   Tsys = 150 + 290*(10^(1.2/10) - 1) = 242.3 K  ->  35.3 - 23.84 = 11.46 -> 11.5 dB/K
cfg.radio.GT_dBK             = 11.5;     % dB/K G/T terminal VSAT Ku (MAESTRO)  [cadena de ruido TR 38.821 Tab. 6.1.1.1-3]
cfg.radio.CNreq_dB           = 5;        % dB   C/N requerido (umbral de enlace)      % <-- DECISION
cfg.radio.cellOffset_deg     = 0;        % deg  offset angular usuario-centro de celda
cfg.radio.p618_availability  = 0.1;      % %    excedencia P.618 lluvia (0.1% = 99.9% disponibilidad)

% ----- Terminal RX en tierra (necesario para discriminar interferentes) -----
% EN ESTE PERFIL G/T es el parametro MAESTRO; Tsys se DERIVA para no romper el C/N:
%   GT = Grx_max - 10*log10(Tsys)  =>  Tsys = 10^((Grx_max - GT)/10)
% Asi el ruido del modulo de interferencia es coherente con compute_link_budget.
%
% OJO (auditoria, hallazgo G3): el MAESTRO NO es el mismo en los dos perfiles y las
% dos elecciones son correctas, luego Tsys_K NO se puede derivar en el punto de uso
% como se hace con B_Hz:
%   config_default (Ku, estudio)  -> G/T MAESTRO,   Tsys derivado   (aqui)
%   config_calib   (Ka, E0 38.821) -> Tsys MAESTRO (Tant + T0*(F-1), Tabla 6.1.1.1-3),
%                                     G/T derivado
% Lo que SI es comun a ambos es el invariante  GT_dBK == Grx_max_dBi - 10log10(Tsys_K),
% y eso es lo que verifica assert_cfg_coherent (llamado desde compute_link_budget y
% run_one_density): si alguien sobrescribe uno sin recalcular el otro, ABORTA en vez de
% devolver un ruido termico equivocado con pinta de correcto.
% Ganancia REAL del terminal: reflector de 60 cm a 12 GHz.
%   lambda = c/f = 0.024983 m  ->  D/lambda = 0.60/0.024983 = 24.02
%   ITU-R S.1428-1:  G = 20*log10(D/lambda) + 7.7 = 27.61 + 7.7 = 35.31 dBi
%   Coherente con la formula de apertura para eta = 0.60:
%     G = 10*log10(eta*(pi*D/lambda)^2) = 35.33 dBi   (coincide, ver comprobacion
%     que imprime run_ffr_demo / test_terminal_coherence)
cfg.radio.Grx_max_dBi          = 35.3;   % dBi  terminal 60 cm @ 12 GHz  [ITU-R S.1428-1]
cfg.radio.Tsys_K               = 10^((cfg.radio.Grx_max_dBi - cfg.radio.GT_dBK)/10);  % K (DERIVADO de G/T en este perfil, NO hardcodear)
cfg.radio.term_efficiency      = 0.6;    % -    eficiencia de apertura del terminal (deriva D/lambda en S.1428)  % <-- DECISION
cfg.radio.term_beamwidth3dB_deg = 2.0;   % deg  ancho de haz a -3 dB del terminal VSAT (NO USADO: el terminal usa S.1428)  % <-- DECISION
cfg.radio.sidelobe_floor_dB    = -30;    % dB   suelo de lobulos del HAZ del satelite (rel. al pico; S.1528 = futuro)  % <-- DECISION

% ----- Interferencia co-canal -----
cfg.interf.reuse      = 1;            % 1 = reuso pleno (peor caso; la FFR lo sobrescribira)   % <-- DECISION
cfg.interf.satPointing = 'nadir';    % 'nadir' (nominal: haz a su nadir) | 'boresight' (cota superior: haz al usuario)  % <-- DECISION

% ----- Layout MULTIHAZ (celdas Earth-fixed; ver build_beam_layout) -----
% Dos formas EXCLUYENTES de fijar el tamano del cluster:
%   nRings = R  -> anillos hex completos, nBeams = 3*R*(R+1)+1 (2->19, 5->91)
%   nRings = [] -> se usa nBeams tal cual  (es el caso de la CALIBRACION E0)
cfg.beams.nBeams   = 19;              % haces del cluster (centro + 2 anillos)    % <-- DECISION
cfg.beams.nRings   = [];              % [] = usar nBeams (E0 depende de esto)
cfg.beams.nUsers   = 1000;            % usuarios sinteticos en el haz central (solo E0/compute_beam_cir)
cfg.beams.cellFrac = 0.5;             % radio de celda = cellFrac*s (crossover en s/2)
cfg.beams.h_ref_km = [];              % [] = altitud de la constelacion 1 (proyeccion de celdas)

% ----- FFR (reutilizacion fraccionaria de frecuencias) -----
% Motor: ffr_coloring -> ffr_policy -> ffr_allocate -> compute_sinr_ffr -> compute_kpis
cfg.ffr.scheme        = 'ffr';        % 'reuse1' | 'reuseD' | 'ffr'                % <-- DECISION
cfg.ffr.Delta         = 3;            % numero de sub-bandas de borde (colores): 3 o 4  % <-- DECISION
cfg.ffr.alpha         = 0.4;          % fraccion de banda del INTERIOR: B = alpha*B + Delta*(1-alpha)*B/Delta  % <-- DECISION
cfg.ffr.classifier    = 'sinr';       % 'sinr' (SINR_ref reuso-1 >= tau_dB) | 'theta' (offset <= tau_theta_deg)
cfg.ffr.tau_dB        = 5;            % dB   umbral centro/borde del clasificador 'sinr'   % <-- DECISION
cfg.ffr.tau_theta_deg = 0.5;          % deg  umbral centro/borde del clasificador 'theta'  % <-- DECISION
% Como se fija el umbral tau (lo decide ffr_policy):
%   'abs'      -> valor absoluto de tau_dB / tau_theta_deg (arriba).
%   'quantile' -> tau se calcula POR INSTANTE para que la fraccion de CENTRO sea
%                 cfg.ffr.tau_q %. Es robusto: el nivel absoluto de SINR en reuso-1
%                 depende mucho de la geometria (compresion angular del cluster
%                 Earth-fixed a baja elevacion), asi que un tau fijo puede dejar
%                 el 100% de usuarios en un lado y degenerar la FFR.
cfg.ffr.tau_mode      = 'abs';        % 'abs' | 'quantile'                        % <-- DECISION
cfg.ffr.tau_q         = 50;           % %  fraccion objetivo de usuarios de CENTRO (modo 'quantile')

% ----- FFR ADAPTATIVA (H2, experimento E5; regla implementada en ffr_policy) -----
% alpha(t) se interpola LINEALMENTE (con saturacion) entre dos anclas de elevacion:
%     w(t)     = clip((elev(t)-elMin_deg)/(elMax_deg-elMin_deg), 0, 1)
%     alpha(t) = alpha_min + (alpha_max-alpha_min)*w(t)      (adapt_dir = +1)
% OJO CON EL SENTIDO: aqui alpha es la fraccion del INTERIOR, luego MAS alpha =
% MENOS banda al BORDE. Proteger el borde en geometria mala = alpha BAJO a baja
% elevacion, que es lo que hace adapt_dir = +1 (alpha CRECIENTE con la elevacion).
% Es ademas el sentido MEDIDO: en el barrido de alpha fijo el optimo pasa de 0.5
% (ventana completa) a 0.6 (elevacion >= 45 deg). adapt_dir = -1 invierte la regla
% y se conserva solo para poder CONTRASTAR el sentido empiricamente en E5.
cfg.ffr.adaptive      = false;        % true = alpha(t) adaptativa (H2)           % <-- DECISION
cfg.ffr.alpha_min     = 0.3;          % alpha en el ancla de elevacion BAJA (geometria mala: protege borde)  % <-- DECISION
cfg.ffr.alpha_max     = 0.7;          % alpha en el ancla de elevacion ALTA (geometria buena: libera agregado)  % <-- DECISION
cfg.ffr.elMin_deg     = 25;           % deg  ancla inferior (elevacion de mascara)  % <-- DECISION
cfg.ffr.elMax_deg     = 90;           % deg  ancla superior (cenit)                 % <-- DECISION
cfg.ffr.adapt_dir     = +1;           % +1 = alpha CRECE con la elevacion (medido) | -1 = inversa
% PENDIENTE (extension): adaptar tambien el cuantil tau (cfg.ffr.tau_q) segun
% POL.edgeFrac_ref / POL.covFrac. Gancho descrito y comentado en ffr_policy.

% Convencion de planificacion (ver ffr_allocate; es METODOLOGICAMENTE CRITICA):
%   'share' -> reparto equitativo de la sub-banda entre los usuarios de la celda que
%              la usan. Es la convencion justa: la FFR parte la poblacion en dos
%              grupos y su sub-banda la comparten MENOS usuarios.
%   'peak'  -> cada usuario recibe la sub-banda completa (tasa asignable de pico).
%              SESGA la comparacion contra la FFR; se deja como diagnostico.
cfg.ffr.sched         = 'share';      % 'share' | 'peak'                          % <-- DECISION

% ----- COMPUTO (rendimiento; NO afecta a ningun resultado numerico) -----
% Troceado TEMPORAL del pipeline (ver run_one_density). El pico de RAM de un punto
% no son los tensores de geometria sino el interior de compute_interference, donde
% coexisten G.* + FSPL_all + Latm + temporales de interp1 sobre arrays [M x N x Nt]
% (~74 B por elemento). Procesar la ventana en BLOQUES de instantes hace que el pico
% lo fije el TAMANO DE BLOQUE y no Nt completo, sin cambiar ningun resultado: el
% pipeline no acopla instantes distintos (geometria, tau por cuantil, carga n0 y
% alpha son por instante; los KPIs son percentiles sobre la union de muestras).
%   []  o  Inf  -> sin trocear (comportamiento por defecto, ventana entera de golpe)
%   n           -> n instantes por bloque
cfg.compute.timeBlock = [];           % instantes por bloque                      % <-- DECISION

% ----- ACOPLAMIENTO POTENCIA-DIRECTIVIDAD AL BARRER EL ANCHO DE HAZ (eje de E3) -----
% Vive AQUI, en el perfil, y no en run_e3_beamdensity: es una hipotesis FISICA sobre
% el sistema (como reparte el satelite su potencia al cambiar el tamano de celda),
% no un parametro del diseno del experimento. Cualquier runner que barra
% cfg.radio.beamwidth3dB_deg debe declararla, no solo E3.
%
% Al estrechar el haz de theta_ref a theta la directividad sube 20*log10(theta_ref/theta)
% (directividad ~ 1/theta^2). Que ocurre con la POTENCIA por haz:
%   'total_const'    -> potencia TOTAL del satelite constante. El area se tesela con
%                       ~1/theta^2 haces, luego la potencia POR HAZ baja en el mismo
%                       factor y las dos correcciones SE CANCELAN:
%                       EIRPdens_eff = EIRPdens_ref  => C/N INVARIANTE al densificar.
%                       Es lo que permite atribuir TODO lo medido en E3 a la
%                       INTERFERENCIA, que es lo que se queria aislar.
%   'per_beam_const' -> potencia por haz constante: EIRPdens_eff = EIRPdens_ref + dG.
%                       Sensibilidad; contrastado en E3b (max 0.277 dB de SINR_edge_p5,
%                       veredictos identicos: el sistema esta limitado por interferencia).
%
% OJO: el acoplamiento viaja por EIRPdensity_dBWMHz, NO por
% Gmax_dBi, que es informativo y se cancela exacto en compute_link_budget.
cfg.e3.powerMode = 'total_const';     % 'total_const' | 'per_beam_const'          % <-- DECISION

% ----- KPIs -----
cfg.kpi.gamma0_dB = 0;                % dB   umbral de cobertura: P(SINR >= gamma0)  % <-- DECISION

% ----- CRITERIO DE VIABILIDAD (H1) -----
% Formaliza el "umbral de viabilidad" que hasta ahora era una tarea abierta. Lo
% evalua compute_kpis y lo devuelve en KPI.viab.* (metrica, cobertura y veredicto
% de 3 niveles inviable / marginal / viable).
%
% TRAZABILIDAD del criterio (para la defensa):
%   - Percentil 5 de los usuarios de BORDE = definicion normalizada de "cell edge"
%     de ITU-R M.2135-1 (Guidelines IMT-Advanced, "5% user spectral efficiency") y
%     de ITU-R M.2410-0 (requisitos minimos IMT-2020, "5th percentile user spectral
%     efficiency"). 3GPP la hereda en la metodologia de evaluacion NTN de
%     TR 38.821, que es EXACTAMENTE la usada en la calibracion E0 de este TFG.
%   - Suelo de SINR = umbral de demodulacion del MCS mas robusto de 5G NR
%     (QPSK, tasa de codigo mas baja, criterio de 10% BLER; tabla MCS/CQI de
%     3GPP TS 38.214). Por debajo de ese suelo NO hay formato de modulacion y
%     codificacion que cierre el enlace: la capacidad de Shannon podra ser >0 pero
%     el sistema real no transmite -> region INVIABLE.
%
% INDEPENDENCIA DE BANDA (por que el criterio es transportable Ka -> Ku):
%   El criterio se formula sobre la SINR de borde y sobre umbrales de MODULACION,
%   magnitudes independientes de la portadora. Todo lo que depende de la banda
%   (FSPL, atenuacion por lluvia y gases, ganancias de apertura) ya esta DENTRO
%   del calculo de la SINR; una vez obtenida la SINR, el umbral que decide si el
%   enlace cierra lo fija la modulacion/codificacion (QPSK 5G NR), no la frecuencia.
%   Por eso, aunque la calibracion E0 use los parametros Ka de TR 38.821 (unica
%   configuracion NTN con resultados publicados con los que contrastar) y el
%   escenario de estudio sea Ku, el criterio se aplica sin reajuste: E0 valida el
%   MOTOR (COMO se calcula la SINR), no la banda concreta.
cfg.viab.metric         = 'sinr_edge_p5';  % metrica reina: percentil 5 de SINR de BORDE  % <-- DECISION
cfg.viab.gamma_floor_dB = -6.7;   % dB  suelo FISICO: MCS0/CQI1 QPSK 5G NR (10% BLER)     % <-- DECISION
cfg.viab.gamma_th_dB    =  0.0;   % dB  umbral primario/conservador (outage satelite)     % <-- DECISION
cfg.viab.cov_target     =  0.95;  % -   P(SINR_borde >= gamma_th) exigida (cobertura 95%) % <-- DECISION
cfg.viab.se_target_bhz  =  0.12;  % bit/s/Hz  objetivo de SE_p5 Rural (ITU-R M.2410-0 Tabla 1)  % <-- DECISION

end
