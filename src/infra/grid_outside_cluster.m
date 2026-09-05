function [isOut, pct, thMax_deg] = grid_outside_cluster(theta_max_deg, spacing_deg, lattice)
%GRID_OUTSIDE_CLUSTER  Criterio UNICO de "la rejilla sobresale del cluster de haces".
%   isOut                       = grid_outside_cluster(theta_max_deg, spacing_deg)
%   [isOut, pct]                = grid_outside_cluster(..., lattice)
%   [isOut, pct, thMax_deg]     = grid_outside_cluster(...)
%
%   lattice = 'hex' (defecto) | 'rect'. Ver "DEPENDENCIA DE LA TOPOLOGIA" abajo: el
%   circunradio de la celda NO es el mismo en las dos retículas, y usar el hexagonal
%   sobre una cuadrada produce FALSOS POSITIVOS. Omitir el argumento reproduce el
%   comportamiento anterior BIT-IDENTICO (todos los llamadores historicos son 'hex').
%
%   FUENTE UNICA del criterio (mismo principio que mem_peak_model_GB y merge_cfg:
%   una sola implementacion, reutilizada). Lo usan:
%     - ffr_context          -> aviso operativo del MOTOR, en cada simulacion
%     - run_e3b_minelev      -> senal de sesgo por punto del barrido de minElev
%     - run_convergence_nrings -> campo R.gridOut del barrido de tamano de cluster
%   Tener el umbral escrito tres veces es justo el patron "dos reglas para lo mismo"
%   que ya causo los fallos de atm_key_of (merge superficial) y de associate_serving
%   (servidor de otra constelacion): en cuanto una copia se corrige y las otras no,
%   el motor y el diagnostico dejan de contar la misma historia.
%
%   GEOMETRIA DEL UMBRAL
%   --------------------
%   Las celdas teselan el suelo en hexagonos con separacion s = spacing_deg entre
%   centros (regla normativa 3GPP TR 38.821 Sec. 6.1.1, ver build_beam_layout). En
%   un hexagono regular con esa separacion entre centros, el punto MAS ALEJADO de su
%   propio centro es un VERTICE, situado a
%
%       circunradio = s / sqrt(3)
%
%   (la apotema vale s/2 y el circunradio apotema*2/sqrt(3) = s/sqrt(3)). Por tanto
%   un usuario cuyo offset angular a su haz servidor supere s/sqrt(3) esta, por pura
%   geometria, FUERA de la celda que se le ha asignado: la rejilla de usuarios
%   sobresale del cluster de haces, esos usuarios se asocian a un haz muy desapuntado
%   y ven MENOS interferentes co-canal de los que les corresponderian -> su SINR es
%   OPTIMISTA y el punto NO es defendible.
%
%   CORRECCION (auditoria, hallazgo M6). Antes las tres copias comparaban contra s,
%   que es sqrt(3) = 1.73 veces mas laxo y dejaba ~2x de zona muerta. Medido sobre
%   los 34 puntos ya publicados (E3 A/B, E3b, nRings, run_ffr_demo/E5) la separacion
%   con el umbral correcto es limpia y no depende de donde se ponga la raya:
%       puntos SANOS    :  92.5% - 98.4% del circunradio (saturan por debajo, que es
%                          lo que predice la geometria: el peor usuario cae en un
%                          vertice de celda)
%       puntos SESGADOS : 117.7% - 904.4%
%   Con el umbral antiguo los sanos ocupaban solo el 53-57%, y en esa zona muerta se
%   colaron dos FALSOS NEGATIVOS (E3 variante A con HPBW=1.0 deg, y nRings=3): ambos
%   los habia marcado ya la guarda ANALITICA de anillos, luego ninguna conclusion
%   publicada estaba afectada, pero los dos criterios se contradecian.
%
%   La justificacion que acompanaba al umbral antiguo ("mas un margen por la
%   compresion angular a baja elevacion") estaba AL REVES y conviene dejarlo escrito:
%   la compresion angular ENCOGE theta visto desde el satelite, luego hace el chequeo
%   MENOS sensible, no mas. Si hiciera falta margen seria en el sentido contrario.
%
%   CRITERIO COMPLEMENTARIO (no lo cubre esta funcion). Existe ademas una guarda
%   ANALITICA de anillos, calculada solo con el layout y sin simular:
%       guarda = (dist. del anillo externo - radio de la rejilla) / s   [anillos]
%   Son criterios INDEPENDIENTES: la guarda analitica mira el layout ANTES de
%   simular; esta funcion mira lo que de verdad ha ocurrido (el offset maximo
%   medido). Mantener los dos es lo que permitio detectar la discrepancia de M6.
%
%   DEPENDENCIA DE LA TOPOLOGIA (anadido al introducir el modo 'rect' de
%   build_beam_layout para el contraste de OneWeb)
%   ---------------------------------------------------------------------------
%   El circunradio de la celda depende de la RETICULA, no solo de s:
%       hexagonal : el punto mas lejano es un vertice, a  s/sqrt(3)  = 0.5774*s
%       cuadrada  : el punto mas lejano es una esquina, a s*sqrt(2)/2 = 0.7071*s
%   Aplicar el umbral hexagonal a una reticula cuadrada es 1.2247 veces demasiado
%   ESTRICTO y produce FALSOS POSITIVOS. Se observo al ejecutar la variante 4x4 de
%   OneWeb: max offset 9.2854 deg avisaba como "fuera de celda" contra 7.7135 deg
%   (120.38%) cuando el limite correcto de una celda cuadrada es 9.4471 deg (98.29%,
%   VALIDO).
%     [cifras CORREGIDAS contra oneweb_results.mat (OUTR.diag.theta_max_deg =
%      9.285378). Antes figuraban aqui 9.118 / 96.5% / 118.2%, que eran erroneas; la
%      conclusion no cambia: >100% con el umbral hexagonal, <100% con el correcto.]
%   Es el mismo modo de fallo que M6 pero con el signo contrario: alli el
%   umbral era demasiado laxo y colaba puntos sesgados; aqui es demasiado estricto y
%   marcaria puntos sanos. En ambos casos la causa es la misma -- un umbral
%   geometrico aplicado fuera de la geometria para la que se dedujo.
%
%   Entradas:
%     theta_max_deg - offset angular MAXIMO medido de un usuario a su haz servidor
%                     [deg] (= max(ctx.theta_deg(ctx.cov)) o diag.theta_max_deg)
%     spacing_deg   - separacion inter-haz s del layout [deg] (BL.spacing_deg)
%     lattice       - (opc.) 'hex' (defecto) | 'rect'. Normalmente BL.lattice.
%
%   Salidas:
%     isOut      - logico: true si la rejilla sobresale (punto SESGADO)
%     pct        - theta_max / circunradio en %, la magnitud INFORMATIVA: cuanto del
%                  circunradio ocupa el peor usuario. 100% = justo en el vertice.
%     thMax_deg  - el umbral aplicado [deg]
%
%   Aportacion propia (infraestructura de diagnostico; NO interviene en ningun
%   calculo fisico y no cambia ningun KPI).

if nargin < 3 || isempty(lattice), lattice = 'hex'; end
switch lower(lattice)
    case 'hex',  k = 1/sqrt(3);        % vertice del hexagono
    case 'rect', k = sqrt(2)/2;        % esquina del cuadrado
    otherwise
        error('grid_outside_cluster:lattice', ...
            'lattice = ''%s'' no soportado (usa ''hex'' o ''rect'').', lattice);
end

thMax_deg = k * spacing_deg;
isOut     = theta_max_deg > thMax_deg;
pct       = 100 * theta_max_deg / thMax_deg;
end
