function MEM = estimate_sweep_memory(cfg_base, pts, nWorkers, opt)
%ESTIMATE_SWEEP_MEMORY  Dimensiona la memoria de un barrido ANTES de lanzarlo.
%   MEM = estimate_sweep_memory(cfg_base, pts)
%   MEM = estimate_sweep_memory(cfg_base, pts, nWorkers)
%   MEM = estimate_sweep_memory(cfg_base, pts, nWorkers, opt)
%
%   Responde a la unica pregunta que decide si un barrido se puede lanzar en esta
%   maquina: CUANTOS WORKERS CABEN. Cada worker de parfor es un proceso con su
%   PROPIA geometria (M*N*Nt), no se comparte nada, luego el coste total es el pico
%   de un punto MULTIPLICADO por el numero de workers, y el limite real de E3 es la
%   RAM, no la CPU (ver "Memoria y troceado temporal" en el README).
%
%   DOS CORRECCIONES respecto al estimador que vivia dentro de run_sweep_points:
%
%   (a) MODELO REAL DE 74 B/elem (antes 25 B/elem, 3x POR DEBAJO). El pico NO son
%       los tensores de compute_geometry (G.el/az/range + G.vis = 25 B/elem), que
%       es lo unico que contaba el estimador anterior: se alcanza DENTRO de
%       compute_interference, en atm_loss_dB(cfg, G.el), donde sobre el mismo array
%       [M x N x Ntb] coexisten ademas FSPL_all, Latm_clear/rain, valid y los
%       temporales de interp1 -> ~74 B/elem, mas ~40 B/elem de ffr_context sobre
%       [M x nBeams x Ntb]. El modelo vive en mem_peak_model_GB (fuente unica,
%       compartida con run_one_density).
%
%   (b) TROCEADO TEMPORAL (cfg.compute.timeBlock). El pico lo fija el tamano de
%       BLOQUE, no la ventana completa: pico ~ 74*M*N*nBlk con nBlk =
%       min(timeBlock, Nt). El estimador anterior ignoraba el troceado y por tanto
%       no servia para lo unico que el troceado existe: decidir cuantos workers
%       caben. Aqui se resuelve nBlk con la MISMA regla que run_one_density.
%
%   Entradas:
%     cfg_base - configuracion base del barrido (M, Nt, nBeams y timeBlock salen
%                de aqui; solo la densidad N varia por punto).
%     pts      - cell {1 x nP} de structs de punto (o un struct suelto). De cada
%                punto se lee la densidad igual que hace run_one_density:
%                pt.constellations (suma de T) tiene prioridad sobre pt.T, y en su
%                defecto se usa cfg_base.constellations(1).T.
%     nWorkers - workers previstos (defecto 1). Solo se usa para el total y para la
%                recomendacion de timeBlock; no altera el pico por worker.
%     opt      - struct opcional:
%                .verbose      false (defecto) = no imprime nada; true = tabla.
%                .budget_frac  fraccion de la RAM libre que se considera utilizable
%                              (defecto 0.7; el 30% restante es margen para el SO,
%                              MATLAB cliente y la fragmentacion del heap).
%                .budget_GB    presupuesto ABSOLUTO en GB; si se da, ignora la RAM
%                              libre medida (util para dimensionar una maquina
%                              distinta de esta, p.ej. un cluster).
%
%   Salidas (struct MEM):
%     .bytes_per_elem   74 (modelo M*N*Ntb) y .bytes_per_elem_beams 40
%     .M .Nt .nBeams    tamanos del escenario
%     .timeBlock_cfg    cfg.compute.timeBlock tal cual ([] / Inf = sin trocear)
%     .nBlk .nBlocks    instantes por bloque efectivos y numero de bloques
%     .N                [1 x nP] satelites de cada punto
%     .peak_GB          [1 x nP] pico por worker CON el troceado aplicado
%     .peak_full_GB     [1 x nP] pico si se hiciera la ventana de una pieza
%     .peak_worker_GB   peor caso de .peak_GB (el que dimensiona)
%     .freeRAM_GB       RAM fisica libre medida (NaN si el SO no la expone)
%     .budget_GB        presupuesto usado en la decision
%     .nWorkers_req     workers previstos
%     .total_GB         peak_worker_GB * nWorkers_req
%     .nFit             workers que CABEN en el presupuesto (>= 1, NaN sin medida)
%     .fits             true si nWorkers_req <= nFit
%     .timeBlock_reco   timeBlock necesario para que quepan nWorkers_req workers
%                       ([] si con el actual ya caben, NaN si no se pudo decidir)
%
%   COHERENCIA con lo MEDIDO (FASE A, config convergida step_km=4 -> M=317,
%   Nt=61, nBeams=91): T=1584 sin trocear -> 2.18 GB (medido 2.18) y con
%   timeBlock=20 -> 0.71 GB (medido 0.71). El estimador anterior daba 0.72 GB para
%   la ventana ENTERA, es decir el valor correcto del bloque atribuido a algo 3x
%   mayor. Verificado en test_mem_estimator.
%
%   Aportacion propia (infraestructura de experimentos; NO interviene en ningun
%   calculo fisico ni cambia ningun resultado de simulacion).

if nargin < 3 || isempty(nWorkers), nWorkers = 1; end
if nargin < 4 || isempty(opt), opt = struct(); end
verbose     = getf(opt, 'verbose',     false);
budgetFrac  = getf(opt, 'budget_frac', 0.7);
budgetGBopt = getf(opt, 'budget_GB',   []);

if ~iscell(pts), pts = {pts}; end
nP = numel(pts);
nWuse = max(1, round(nWorkers));

%% 1. Tamanos del escenario (no dependen del punto: solo varia la densidad)
users = build_user_grid(cfg_base);
BL    = build_beam_layout(cfg_base);
M     = users.M;
nB    = BL.nBeams;
Nt    = numel(cfg_base.time.t0 : cfg_base.time.dt : cfg_base.time.t0 + cfg_base.time.duration);

% Bloque temporal EFECTIVO, con la misma regla que run_one_density (replica de
% diagnostico, igual que atm_key_of replica la clave de la cache P.618: si la regla
% del motor cambiase habria que cambiarla tambien aqui, por eso se documenta).
[nBlk, nBlocks, tbCfg] = resolve_time_block_diag(cfg_base, Nt);

%% 2. Densidad de cada punto y pico por worker
N = zeros(1, nP);
for i = 1:nP
    Ni = cfg_base.constellations(1).T;
    if isfield(pts{i},'T') && ~isempty(pts{i}.T), Ni = pts{i}.T; end
    if isfield(pts{i},'constellations') && ~isempty(pts{i}.constellations)
        Ni = sum([pts{i}.constellations.T]);      % saturacion INTER-constelacion
    end
    N(i) = Ni;
end

peak_GB      = mem_peak_model_GB(M, N, nBlk, nB);   % con troceado (el que manda)
peak_full_GB = mem_peak_model_GB(M, N, Nt,   nB);   % de una pieza (referencia)
peakW        = max(peak_GB);

%% 3. Presupuesto de RAM y workers que caben
freeRAM_GB = free_ram_GB();
if ~isempty(budgetGBopt)
    budget_GB = budgetGBopt;
else
    budget_GB = budgetFrac * freeRAM_GB;            % NaN si no hay medida
end

if isfinite(budget_GB) && budget_GB > 0
    nFit = max(1, floor(budget_GB / max(peakW, eps)));
else
    nFit = NaN;
end
total_GB = peakW * nWuse;
fits     = ~isfinite(nFit) || nWuse <= nFit;

% timeBlock necesario para que quepan los workers pedidos. El pico es LINEAL en
% nBlk (pico = M*nBlk*(N*74 + nB*40)), asi que se despeja directamente.
timeBlock_reco = [];
if ~fits && isfinite(budget_GB)
    perInst = M * (max(N)*74 + nB*40) / 2^30;       % GB por instante de bloque
    tbr = floor((budget_GB / nWuse) / max(perInst, eps));
    if tbr < 1
        timeBlock_reco = NaN;                       % ni 1 instante cabe
    else
        timeBlock_reco = min(tbr, Nt);
    end
end

%% 4. Salida
MEM = struct( ...
    'bytes_per_elem',       74, ...
    'bytes_per_elem_beams', 40, ...
    'M', M, 'Nt', Nt, 'nBeams', nB, ...
    'timeBlock_cfg',  tbCfg, ...
    'nBlk', nBlk, 'nBlocks', nBlocks, ...
    'N', N, ...
    'peak_GB',        peak_GB, ...
    'peak_full_GB',   peak_full_GB, ...
    'peak_worker_GB', peakW, ...
    'freeRAM_GB',     freeRAM_GB, ...
    'budget_GB',      budget_GB, ...
    'budget_frac',    budgetFrac, ...
    'nWorkers_req',   nWuse, ...
    'total_GB',       total_GB, ...
    'nFit',           nFit, ...
    'fits',           fits, ...
    'timeBlock_reco', timeBlock_reco);

if verbose, print_mem_table(MEM); end
end

% -------------------------------------------------------------------------
function print_mem_table(MEM)
%PRINT_MEM_TABLE  Tabla de memoria del barrido.
fprintf('[memoria] modelo %d B/elem (M*N*Nt) + %d B/elem (M*nBeams*Nt) -- FASE A, cota superior\n', ...
    MEM.bytes_per_elem, MEM.bytes_per_elem_beams);
if MEM.nBlocks > 1
    fprintf('          M=%d usuarios | Nt=%d instantes | nBeams=%d | troceado: timeBlock=%d -> %d bloques\n', ...
        MEM.M, MEM.Nt, MEM.nBeams, MEM.nBlk, MEM.nBlocks);
else
    fprintf('          M=%d usuarios | Nt=%d instantes | nBeams=%d | SIN trocear (cfg.compute.timeBlock vacio)\n', ...
        MEM.M, MEM.Nt, MEM.nBeams);
end

% Con muchos puntos solo interesan los mas caros: el pico es monotono en N.
nP = numel(MEM.N);
[~, ord] = sort(MEM.N, 'descend');
show = ord(1:min(nP, 8));
show = sort(show);
fprintf('          %-6s %8s %16s %16s\n', 'punto', 'N sat', 'pico/worker[GB]', 'sin trocear[GB]');
for k = 1:numel(show)
    i = show(k);
    fprintf('          %-6s %8d %16.2f %16.2f\n', sprintf('p%02d', i), MEM.N(i), ...
        MEM.peak_GB(i), MEM.peak_full_GB(i));
end
if numel(show) < nP
    fprintf('          (%d puntos mas, todos por debajo del peor caso)\n', nP - numel(show));
end

if isfinite(MEM.freeRAM_GB)
    fprintf('          RAM libre %.2f GB | presupuesto %.0f%% = %.2f GB | pico/worker %.2f GB -> caben %d worker(s)\n', ...
        MEM.freeRAM_GB, 100*MEM.budget_frac, MEM.budget_GB, MEM.peak_worker_GB, MEM.nFit);
else
    fprintf('          RAM libre no disponible en esta plataforma | pico/worker %.2f GB (sin veredicto)\n', ...
        MEM.peak_worker_GB);
end

if MEM.fits
    if isfinite(MEM.nFit)
        fprintf('          %d worker(s) previstos x %.2f GB = %.2f GB -> CABE\n', ...
            MEM.nWorkers_req, MEM.peak_worker_GB, MEM.total_GB);
    end
else
    fprintf(['[AVISO] %d worker(s) previstos x %.2f GB = %.2f GB > presupuesto %.2f GB.\n' ...
             '        Cada worker mantiene su PROPIA geometria, no se comparte.\n'], ...
        MEM.nWorkers_req, MEM.peak_worker_GB, MEM.total_GB, MEM.budget_GB);
    if isempty(MEM.timeBlock_reco)
        % nada que recomendar
    elseif isnan(MEM.timeBlock_reco)
        fprintf(['        Ni con timeBlock=1 caben: reduce nWorkers (caben %d), la rejilla de\n' ...
                 '        usuarios (cfg.ground.step_km) o la densidad del punto mas caro.\n'], MEM.nFit);
    else
        fprintf(['        Opciones: usar %d worker(s), o fijar cfg.compute.timeBlock = %d\n' ...
                 '        (pico/worker -> %.2f GB), que es EXACTO y no cambia ningun resultado.\n'], ...
            MEM.nFit, MEM.timeBlock_reco, ...
            mem_peak_model_GB(MEM.M, max(MEM.N), MEM.timeBlock_reco, MEM.nBeams));
    end
end
end

% -------------------------------------------------------------------------
function [nBlk, nBlocks, tbCfg] = resolve_time_block_diag(cfg, Nt)
%RESOLVE_TIME_BLOCK_DIAG  Replica de diagnostico de resolve_time_block.
%   MISMA regla que la subfuncion homonima de run_one_density: [] / Inf / <=0 /
%   >= Nt -> un solo bloque. Se replica (no se comparte) por el mismo motivo que
%   atm_key_of replica la clave de la cache P.618: es diagnostico previo al bucle y
%   no debe poder alterar el camino de ejecucion del motor. Si la regla del motor
%   cambiara, hay que cambiarla tambien aqui.
tbCfg = [];
tb = Inf;
if isfield(cfg,'compute') && isfield(cfg.compute,'timeBlock') && ~isempty(cfg.compute.timeBlock)
    tb = cfg.compute.timeBlock;
    tbCfg = tb;
end
if ~isfinite(tb) || tb <= 0 || tb >= Nt
    nBlk = Nt;
else
    nBlk = max(1, floor(tb));
end
nBlocks = ceil(Nt / nBlk);
end

% -------------------------------------------------------------------------
function gb = free_ram_GB()
%FREE_RAM_GB  RAM fisica libre [GB], o NaN si el sistema no la expone.
%   `memory` solo existe en Windows; en otras plataformas se devuelve NaN y el
%   veredicto de workers simplemente no se emite (nunca bloquea el barrido).
gb = NaN;
try
    if ispc
        [~, sys] = memory;
        gb = sys.PhysicalMemory.Available / 2^30;
    end
catch
    gb = NaN;
end
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
