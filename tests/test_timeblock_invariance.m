%% TEST_TIMEBLOCK_INVARIANCE  Trocear el tiempo NO debe cambiar ningun resultado.
%  Verificacion de que cfg.compute.timeBlock es una optimizacion de MEMORIA y no un
%  cambio de modelo. Ejecuta EL MISMO punto (run_one_density) con:
%
%     (a) timeBlock = Inf  -> un solo bloque (comportamiento por defecto)
%     (b) timeBlock = 20   -> bloques de 20 instantes
%     (c) timeBlock = 10   -> bloques de 10 instantes
%
%  y compara TODOS los structs de resultados campo a campo: KPIs .all/.center/.edge
%  (medias, percentiles p5/p50/p95, R_agg, Jain, Pcov, SE) MAS el bloque .viab
%  (SINR_edge_p5, SE_edge_p5, cobertura, margenes y veredicto) MAS las CDF
%  (keep_cdf = true), el diagnostico y la definicion de casos.
%
%  Criterio de aceptacion: max|diferencia| == 0 (o < 1e-12 si alguna reduccion
%  cambiase de orden). Cualquier diferencia mayor PARA el test e informa del campo.
%
%  Se prueban los DOS esquemas en el mismo punto (reuse1 y ffr), porque comparten
%  la clasificacion centro/borde de referencia y el troceado tiene que preservar
%  tambien esa dependencia entre esquemas, no solo los numeros de cada uno.
%
%  El Nt del escenario se elige para que los bloques NO dividan exactamente la
%  ventana: con Nt = 41, timeBlock = 20 da bloques [20 20 1] y timeBlock = 10 da
%  [10 10 10 10 1]. El ultimo bloque de UN solo instante es el caso limite donde
%  fallaria cualquier combinacion mal hecha.
%
%  Uso: ejecutar en la raiz del proyecto.

clear; clc;

fprintf('===== TEST DE INVARIANCIA  troceado temporal (cfg.compute.timeBlock) =====\n');

%% 1. Escenario reducido (perfil Ku de estudio)
cfg = config_default();
cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 30;
cfg.ground.step_km   = 4;
cfg.time.dt          = 90;          % s
cfg.time.duration    = 3600;        % s -> Nt = 41 (no divisible por 20 ni por 10)
cfg.ffr.tau_mode     = 'quantile';
cfg.ffr.tau_q        = 50;
cfg.ffr.sched        = 'share';
cfg.beams.nRings     = 3;           % 37 haces

Nt = numel(cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration);

%% 2. Un punto con los DOS esquemas
sch = { struct('scheme','reuse1','Delta',1,'alpha',NaN,'adaptive',false,'name','reuse1'), ...
        struct('scheme','ffr',   'Delta',3,'alpha',0.4,'adaptive',false,'name','ffr(D=3,a=0.4)') };

pt = struct('label','tb', 'T',300, 'P',12, 'cases',{sch}, ...
            'elMin_valid',45, 'keep_cdf',true, 'seed',42, 'idx',1);

fprintf('Escenario: T=%d, Nt=%d instantes, %d usuarios de rejilla, 2 esquemas\n', ...
    pt.T, Nt, numel(build_user_grid(cfg).lat));

%% 3. Precalentar la cache P.618 (coste fijo, fuera de los cronometros)
fprintf('\n[warmup] construyendo la tabla ITU-R P.618...\n');
tw = tic;  atm_loss_dB(cfg, 45);  fprintf('[warmup] listo (%.1f s)\n', toc(tw));

%% 4. Ejecutar con distintos tamanos de bloque
blocks = {Inf, 20, 10};
R = cell(1, numel(blocks));
fprintf('\n%-12s %8s %9s %11s %12s %10s\n', ...
    'timeBlock','bloques','tiempo[s]','pico[GB]','sin trocear','muestras');
fprintf('%s\n', repmat('-',1,66));
for i = 1:numel(blocks)
    cfgB = cfg;
    cfgB.compute.timeBlock = blocks{i};
    R{i} = run_one_density(cfgB, pt);
    fprintf('%-12s %8d %9.1f %11.3f %12.3f %8.1f MB\n', ...
        num2str(blocks{i}), R{i}.mem.nBlocks, R{i}.runinfo.elapsed_s, ...
        R{i}.mem.peak_model_GB, R{i}.mem.peak_full_GB, R{i}.mem.samples_MB);
end
fprintf('%s\n', repmat('-',1,66));

%% 5. COMPARACION contra el caso sin trocear
fprintf('\n===================== COMPARACION =====================\n');
% OUT.runinfo guarda el tiempo de calculo (no determinista) y OUT.mem describe el
% troceado APLICADO, que por definicion cambia entre variantes: ambos se excluyen.
% Todo lo demas -- KPIs, viabilidad, CDF, diagnostico, casos -- debe coincidir.
strip = @(r) rmfield(r, {'runinfo','mem'});
ref   = strip(R{1});

worst = 0;  worstWhere = '-';  allEq = true;
for i = 2:numel(blocks)
    cur = strip(R{i});
    eq  = isequaln(ref, cur);
    [md, where] = deep_maxdiff(ref, cur, 'OUT');
    allEq = allEq && eq;
    if md > worst, worst = md;  worstWhere = where;  end
    fprintf('timeBlock = %-5s (%d bloques): isequaln = %-3s | max|dif| = %.3g   (%s)\n', ...
        num2str(blocks{i}), R{i}.mem.nBlocks, ternary(eq,'SI','NO'), md, where);
end

tol = 1e-12;
if worst > tol || ~isfinite(worst)
    fprintf('\n*** FALLO DE INVARIANCIA ***\n');
    fprintf('El campo %s difiere en %.6g (> %.0e).\n', worstWhere, worst, tol);
    error('test_timeblock_invariance:mismatch', ...
        ['El troceado temporal NO es exacto (campo %s, dif %.6g). PARAR: no usar ' ...
         'cfg.compute.timeBlock en E3 hasta corregirlo.'], worstWhere, worst);
end

%% 6. Valores de referencia (que se ha comparado, no solo que coincide)
fprintf('\nKPIs obtenidos (identicos en las %d variantes):\n', numel(blocks));
fprintf('%-16s %11s %11s %11s %11s %11s\n', ...
    'esquema','SINR_p5_B','SE_p5_B','R_p5_B[Mb]','R_agg[Gb]','veredicto');
fprintf('%s\n', repmat('-',1,78));
for c = 1:numel(sch)
    K = R{1}.K{c};
    fprintf('%-16s %11.4f %11.4f %11.3f %11.4f %11s\n', ...
        R{1}.cases{c}.name, K.viab.SINR_edge_p5, K.viab.SE_edge_p5, ...
        K.edge.R_p5_Mbps, K.all.R_agg_Mbps/1e3, upper(K.viab.verdict));
end
fprintf('%s\n', repmat('-',1,78));

%% 7. Efecto en memoria para la configuracion de E3
%  Mismo modelo que run_one_density/mem_peak_GB (~74 B por elemento M*N*Nt de la
%  geometria + ~40 B por elemento M*nBeams*Nt de ffr_context). Se replica aqui solo
%  para poder informar de la config de E3 sin ejecutarla; la fuente de verdad es
%  run_one_density.
memGB = @(M,N,Ntb,nB) M*Ntb*(N*74 + nB*40) / 2^30;
M_e3 = 317;  N_e3 = 1584;  Nt_e3 = 61;  nB_e3 = 91;    % config convergida (FASE A)
fprintf('\n=========== EFECTO EN E3 (step_km=4 -> %d usuarios, T=%d, Nt=%d, %d haces) ===========\n', ...
    M_e3, N_e3, Nt_e3, nB_e3);
fprintf('%-22s %10s %10s %12s\n','timeBlock','bloques','pico[GB]','reduccion');
fprintf('%s\n', repmat('-',1,58));
full = memGB(M_e3, N_e3, Nt_e3, nB_e3);
for tb = [Nt_e3 30 20 10 5]
    fprintf('%-22s %10d %10.2f %11.2fx\n', ...
        ternary_str(tb >= Nt_e3, 'Inf (sin trocear)', num2str(tb)), ...
        ceil(Nt_e3/tb), memGB(M_e3,N_e3,tb,nB_e3), full/memGB(M_e3,N_e3,tb,nB_e3));
end
fprintf('%s\n', repmat('-',1,58));

%% 8. Veredicto
fprintf('\n===================== VEREDICTO =====================\n');
fprintf('INVARIANCIA: OK  (isequaln = %s en todas las variantes, max|dif| = %.3g <= %.0e)\n', ...
    ternary(allEq,'SI','NO'), worst, tol);
fprintf('Comparados: KPIs .all/.center/.edge (medias, p5/p50/p95, R_agg, Jain, Pcov, SE),\n');
fprintf('            bloque .viab (SINR_edge_p5, cobertura, margenes, veredicto),\n');
fprintf('            CDF completas (keep_cdf = true), .diag y .cases, en 2 esquemas.\n');
fprintf('=====================================================\n');

% =========================================================================
function [md, where] = deep_maxdiff(a, b, path)
%DEEP_MAXDIFF  Maxima diferencia absoluta recorriendo structs/cells recursivamente.
%   Devuelve Inf (y la ruta) si las estructuras no son comparables: distinto tipo,
%   distinto tamano, distintos campos, o NaN en uno y no en el otro.
md = 0;  where = path;

if ~strcmp(class(a), class(b))
    md = Inf;  where = sprintf('%s [clase %s vs %s]', path, class(a), class(b));  return;
end

if isstruct(a)
    if ~isequal(size(a), size(b))
        md = Inf;  where = sprintf('%s [tamano de struct]', path);  return;
    end
    fa = sort(fieldnames(a));  fb = sort(fieldnames(b));
    if ~isequal(fa, fb)
        md = Inf;  where = sprintf('%s [campos distintos]', path);  return;
    end
    for e = 1:numel(a)
        for i = 1:numel(fa)
            if numel(a) > 1
                p = sprintf('%s(%d).%s', path, e, fa{i});
            else
                p = sprintf('%s.%s', path, fa{i});
            end
            [d, w] = deep_maxdiff(a(e).(fa{i}), b(e).(fa{i}), p);
            if d > md, md = d;  where = w;  end
        end
    end
    return;
end

if iscell(a)
    if ~isequal(size(a), size(b))
        md = Inf;  where = sprintf('%s [tamano de cell]', path);  return;
    end
    for i = 1:numel(a)
        [d, w] = deep_maxdiff(a{i}, b{i}, sprintf('%s{%d}', path, i));
        if d > md, md = d;  where = w;  end
    end
    return;
end

if ischar(a) || isstring(a)
    if ~isequal(a, b), md = Inf;  where = sprintf('%s [texto]', path);  end
    return;
end

if isnumeric(a) || islogical(a)
    if ~isequal(size(a), size(b))
        md = Inf;  where = sprintf('%s [tamano %s vs %s]', path, ...
            mat2str(size(a)), mat2str(size(b)));  return;
    end
    if isempty(a), return; end
    x = double(a);  y = double(b);
    d = abs(x - y);
    d(isnan(x) & isnan(y)) = 0;         % NaN vs NaN = igual (isequaln)
    d(xor(isnan(x), isnan(y))) = Inf;   % NaN en uno solo = incomparable
    md = max(d(:));
    return;
end

if ~isequaln(a, b), md = Inf;  where = sprintf('%s [tipo %s]', path, class(a)); end
end

% =========================================================================
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function s = ternary_str(c, a, b)
if c, s = a; else, s = b; end
end
