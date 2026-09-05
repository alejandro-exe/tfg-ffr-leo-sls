function [results, INFO] = run_sweep_points(cfg_base, pts, nWorkers, opt)
%RUN_SWEEP_POINTS  Recorre los puntos de un barrido (secuencial o en paralelo).
%   [results, INFO] = run_sweep_points(cfg_base, pts)
%   [results, INFO] = run_sweep_points(cfg_base, pts, nWorkers)
%   [results, INFO] = run_sweep_points(cfg_base, pts, nWorkers, opt)
%
%   PIEZA DE MOTOR (no runner): el bucle del barrido vive aqui una sola vez y lo
%   reutilizan E3, el Monte Carlo y el test de invariancia, igual que la fisica
%   vive en el motor de la FFR y no en los runners.
%
%   Cada punto es INDEPENDIENTE (run_one_density es una funcion pura: copia local
%   de cfg, sin estado global, sin variables compartidas), luego el bucle es
%   paralelizable sin ninguna reduccion ni seccion critica:
%
%       parfor (i = 1:nP, W)
%           results{i} = run_one_density(cfg_base, pts{i});
%       end
%
%   El cuerpo escribe UNICAMENTE en results{i} (variable de salida troceada por el
%   indice del bucle) y lee cfg_base (variable de difusion) y pts{i} (variable de
%   entrada troceada). No hay ninguna otra escritura fuera del indice.
%
%   nWorkers (clave del diseno):
%     0            -> ejecucion SECUENCIAL. Se usa la forma `parfor (i=1:n, 0)`,
%                     que MATLAB ejecuta en el cliente como un bucle for normal.
%                     Asi el modo secuencial y el paralelo recorren EXACTAMENTE el
%                     mismo camino de codigo y la comparacion de invariancia
%                     (test_parfor_invariance) mide solo el efecto del paralelismo.
%     W > 0        -> hasta W workers (abre el pool si hace falta).
%     [] / omitido -> automatico: todos los workers disponibles, o 0 si no hay
%                     Parallel Computing Toolbox.
%
%   ROBUSTEZ SIN TOOLBOX / SIN WORKERS: `parfor` es sintaxis de MATLAB base. Sin
%   Parallel Computing Toolbox, o con un pool de 0 workers, MATLAB lo ejecuta
%   como un `for` secuencial en el cliente sin error y sin cambiar el resultado.
%   Esta funcion ademas lo detecta antes y degrada explicitamente a nWorkers = 0,
%   dejandolo escrito en INFO.mode en lugar de fallar.
%
%   opt (struct, todos opcionales):
%     .warmup   true (defecto si hay pool) = precalentar la tabla P.618 en cada
%               worker ANTES de cronometrar. Ver la nota de la cache mas abajo.
%     .verbose  true (defecto) = imprimir cabecera y resumen.
%     .mem_budget_frac  fraccion de la RAM libre considerada utilizable (0.7).
%     .mem_budget_GB    presupuesto ABSOLUTO de RAM [GB]; ignora la RAM medida.
%
%   MEMORIA (el limite real de E3, no la CPU): cada worker mantiene su PROPIA
%   geometria M*N*Nt, luego el coste es el pico de UN punto por el numero de
%   workers. La estimacion la hace estimate_sweep_memory con el modelo de la FASE A
%   (74 B/elem) y teniendo en cuenta cfg.compute.timeBlock, que es lo que fija el
%   pico: con troceado el pico lo marca el BLOQUE, no la ventana. Se imprime la
%   tabla (pico/worker, RAM libre, workers que caben) y, si no caben los pedidos,
%   el timeBlock que haria falta. Queda en INFO.mem.
%
%   Salidas:
%     results - cell {1 x nP} con la salida de run_one_density de cada punto,
%               en el ORDEN DE pts (parfor no altera el orden de los resultados,
%               solo el orden de ejecucion).
%     INFO    - struct con el modo, workers, tiempos y diagnostico de la cache.
%
%   ---------------------------------------------------------------------------
%   CACHE PERSISTENT DE atm_loss_dB BAJO PARFOR (comportamiento documentado)
%   ---------------------------------------------------------------------------
%   La tabla ITU-R P.618 de atm_loss_dB se cachea en variables `persistent` con la
%   clave  freq_GHz | lat | lon | minElev | p618_availability.  Cada worker es un
%   proceso MATLAB con su PROPIO espacio de persistents:
%     - construye la tabla la primera vez que se le asigna un punto (~90 s), y
%     - la reutiliza en todos los demas puntos que le toquen.
%   Coste total: una construccion POR WORKER (no por punto), en paralelo. En el
%   modo secuencial: una sola construccion en el cliente.
%   Ninguno de los 5 campos de la clave depende de la densidad de satelites, asi
%   que un barrido de saturacion entero comparte UNA tabla por worker. Esta
%   funcion lo VERIFICA: calcula la clave efectiva de cada punto (misma formula
%   que atm_loss_dB, aqui solo con fines de diagnostico) y avisa si un barrido
%   mezcla configuraciones atmosfericas distintas, que es el unico caso en que la
%   tabla se reconstruiria varias veces por worker.
%
%   Aportacion propia (infraestructura de experimentos; no altera la fisica).

if nargin < 3, nWorkers = []; end
if nargin < 4 || isempty(opt), opt = struct(); end
verbose = getf(opt, 'verbose', true);

%% 1. Guardas sobre la lista de puntos
if ~iscell(pts) || isempty(pts)
    error('run_sweep_points:pts', 'pts debe ser un cell array no vacio de structs de punto.');
end
nP = numel(pts);
for i = 1:nP
    if ~isstruct(pts{i})
        error('run_sweep_points:ptType', 'pts{%d} no es un struct.', i);
    end
    % El indice del punto se fija AQUI, antes del bucle, para que el cuerpo del
    % parfor sea puro: cada punto lleva su propio substream de RNG y su resultado
    % no depende del orden en que los workers lo ejecuten.
    if ~isfield(pts{i},'idx') || isempty(pts{i}.idx), pts{i}.idx = i; end
    if ~isfield(pts{i},'label') || isempty(pts{i}.label)
        pts{i}.label = sprintf('p%02d', i);
    end
end

%% 2. Workers disponibles y degradacion segura
[nEff, pool, hasPCT] = resolve_workers(nWorkers);
if nEff > 0, modeStr = sprintf('parfor (%d workers)', nEff); else, modeStr = 'for (secuencial)'; end

%% 3. Diagnostico: claves de la cache P.618 y memoria estimada por worker
keys = cell(1,nP);
for i = 1:nP
    keys{i} = atm_key_of(cfg_base, pts{i});
end
uk = unique(keys);

% Estimacion de memoria: modelo REAL de 74 B/elem y troceado temporal incluido
% (ver estimate_sweep_memory). El presupuesto de RAM se puede fijar a mano con
% opt.mem_budget_GB / opt.mem_budget_frac para dimensionar otra maquina.
memOpt = struct('verbose', verbose, 'budget_frac', getf(opt,'mem_budget_frac',0.7));
if isfield(opt,'mem_budget_GB') && ~isempty(opt.mem_budget_GB)
    memOpt.budget_GB = opt.mem_budget_GB;
end

if verbose
    fprintf('\n===== BARRIDO: %d puntos | %s =====\n', nP, modeStr);
    fprintf('[cache P.618] %d configuracion(es) atmosferica(s) distinta(s) en el barrido', numel(uk));
    if numel(uk) == 1
        fprintf(' -> cada worker construye la tabla UNA vez y la reutiliza.\n');
    else
        fprintf(' -> la tabla se reconstruira al cambiar de configuracion.\n');
    end
end

MEM = estimate_sweep_memory(cfg_base, pts, max(nEff,1), memOpt);
memPerPoint_GB = MEM.peak_GB;

%% 4. Precalentamiento de la cache P.618 en los workers (NO afecta a los numeros)
%  Solo sirve para que el cronometro mida el barrido y no la construccion de la
%  tabla, que es un coste fijo de arranque de cada worker.
doWarm = getf(opt, 'warmup', nEff > 0);
INFO.warmup_s = 0;
if doWarm && nEff > 0 && hasPCT
    tw = tic;
    try
        cfgw = cfg_base;
        spmd
            atm_loss_dB(cfgw, 45);      % construye/valida la tabla en CADA worker
        end
        INFO.warmup_s = toc(tw);
        if verbose
            fprintf('[warmup] tabla P.618 construida en los %d workers (%.1f s, fuera del cronometro)\n', ...
                nEff, INFO.warmup_s);
        end
    catch ME
        warning('run_sweep_points:warmup', ...
            'No se pudo precalentar la cache P.618 en los workers (%s). Sigue siendo correcto, solo mas lento.', ...
            ME.message);
    end
end

%% 5. EL BUCLE (unico camino de codigo para secuencial y paralelo)
results = cell(1, nP);
tSweep  = tic;

parfor (i = 1:nP, nEff)
    results{i} = run_one_density(cfg_base, pts{i});
end

INFO.elapsed_s = toc(tSweep);

%% 6. Resumen
INFO.mode          = modeStr;
INFO.nPoints       = nP;
INFO.nWorkers_req  = nWorkers;
INFO.nWorkers_eff  = nEff;
INFO.poolSize      = 0;
if ~isempty(pool), INFO.poolSize = pool.NumWorkers; end
INFO.hasPCT        = hasPCT;
INFO.atm_keys      = uk;
INFO.nAtmKeys      = numel(uk);
INFO.mem_per_point_GB = memPerPoint_GB;   % pico por worker CON el troceado aplicado
INFO.mem              = MEM;              % estimacion completa (ver estimate_sweep_memory)

cpu = 0;
for i = 1:nP, cpu = cpu + results{i}.runinfo.elapsed_s; end
INFO.cpu_s = cpu;

if verbose
    fprintf('[barrido] %s: %.1f s de reloj | %.1f s de CPU acumulada (suma de puntos)\n', ...
        modeStr, INFO.elapsed_s, INFO.cpu_s);
end
end

% -------------------------------------------------------------------------
function [nEff, pool, hasPCT] = resolve_workers(nWorkers)
%RESOLVE_WORKERS  Workers realmente utilizables, degradando a 0 sin fallar.
pool = [];
try
    hasPCT = ~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox');
catch
    hasPCT = false;
end

if ~hasPCT
    nEff = 0;                                   % parfor -> for (MATLAB base)
    return;
end

try
    pool = gcp('nocreate');
    if isempty(nWorkers)
        % Automatico: usa el pool existente; si no hay, abre el de por defecto.
        if isempty(pool), pool = parpool('local'); end
        nEff = pool.NumWorkers;
    elseif nWorkers <= 0
        nEff = 0;                               % secuencial EXPLICITO
    else
        if isempty(pool)
            pool = parpool('local', nWorkers);
        elseif pool.NumWorkers < nWorkers
            nWorkers = pool.NumWorkers;         % no se puede pedir mas del pool
        end
        nEff = min(nWorkers, pool.NumWorkers);
    end
catch ME
    warning('run_sweep_points:noPool', ...
        'No hay pool de parfor disponible (%s): el barrido se ejecuta SECUENCIALMENTE.', ME.message);
    pool = [];
    nEff = 0;
end
end

% -------------------------------------------------------------------------
function k = atm_key_of(cfg_base, pt)
%ATM_KEY_OF  Clave de la cache de atm_loss_dB para la cfg efectiva de un punto.
%   DIAGNOSTICO unicamente (no interviene en ningun calculo): replica la formula
%   de clave de atm_loss_dB para poder contar cuantas tablas P.618 distintas
%   necesita el barrido. Si se anadiera un campo a la clave real habria que
%   anadirlo aqui; por eso se comprueba contra los mismos 5 campos.
cfg = cfg_base;
if isfield(pt,'cfg_over') && isstruct(pt.cfg_over)
    % MISMA regla de merge que aplica run_one_density al simular el punto (fichero
    % compartido merge_cfg). Antes se hacia aqui un merge SUPERFICIAL que
    % reemplazaba structs enteros: con un override anidado como
    % .radio.beamwidth3dB_deg se perdian freq_GHz/lat/lon y la clave fallaba.
    cfg = merge_cfg(cfg, pt.cfg_over);
end
k = sprintf('%.10g|%.10g|%.10g|%.10g|%.10g', ...
    cfg.radio.freq_GHz, cfg.ground.point(1), cfg.ground.point(2), ...
    cfg.geom.minElev,   cfg.radio.p618_availability);
end

% -------------------------------------------------------------------------
function v = getf(s, name, defaultVal)
%GETF  Campo opcional de un struct con valor por defecto.
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultVal;
end
end
