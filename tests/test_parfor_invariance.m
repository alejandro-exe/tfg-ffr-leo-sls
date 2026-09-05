%% TEST_PARFOR_INVARIANCE  El barrido en paralelo debe dar resultados IDENTICOS.
%  Verificacion de que paralelizar el barrido con `parfor` es una optimizacion de
%  RENDIMIENTO y no un cambio de modelo. Ejecuta el MISMO barrido pequeno de dos
%  formas y compara los structs de resultados campo a campo:
%
%     (a) SECUENCIAL   run_sweep_points(cfg, pts, 0)   -> parfor con 0 workers,
%                      que MATLAB ejecuta en el cliente como un for normal.
%     (b) PARALELO     run_sweep_points(cfg, pts, [])  -> todos los workers.
%
%  Criterio de aceptacion:
%     - isequaln sobre los resultados (excluyendo OUT.runinfo, que guarda a
%       proposito los metadatos NO deterministas: tiempo de calculo del punto), y
%     - max|diferencia| == 0, o < 1e-12 si alguna suma se hiciera en otro orden.
%  Si aparece CUALQUIER diferencia > 1e-12 el test PARA e informa del campo exacto.
%
%  Escenario REDUCIDO a proposito (rejilla y ventana pequenas, 37 haces): lo que
%  se valida es la INFRAESTRUCTURA del barrido, no la fisica, que ya esta validada
%  en E0 / run_ffr_demo / run_e5_adaptive. Densidades T = [66, 300, 600] con un
%  solo esquema, tal y como pide el diseno del experimento E3.
%
%  Uso: ejecutar en la raiz del proyecto.

clear; clc;

fprintf('===== TEST DE INVARIANCIA  secuencial vs parfor =====\n');

%% 1. Escenario reducido (perfil Ku de estudio)
cfg = config_default();

cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 30;        % rejilla local pequena
cfg.ground.step_km   = 4;

cfg.time.dt       = 600;          % s  (instantes bien separados: mas diversidad
cfg.time.duration = 3600;         % s   de elevacion con muy pocas muestras)

cfg.ffr.tau_mode = 'quantile';    % igual que run_ffr_demo
cfg.ffr.tau_q    = 50;
cfg.ffr.sched    = 'share';

cfg.beams.nRings = 3;             % 37 haces (E3 usara 5 anillos = 91)

%% 2. Puntos del barrido: 3 densidades, UN esquema
Tlist  = [66 300 600];            % T multiplo de P = 6
scheme = struct('scheme','ffr','Delta',3,'alpha',0.4,'adaptive',false,'name','FFR (D=3, a=0.4)');

pts = cell(1, numel(Tlist));
for i = 1:numel(Tlist)
    pts{i} = struct( ...
        'label',       sprintf('T=%d', Tlist(i)), ...
        'T',           Tlist(i), ...
        'P',           6, ...
        'cases',       {{scheme}}, ...
        'elMin_valid', 45, ...
        'keep_cdf',    true, ...     % comparacion MAS estricta: tambien las CDF
        'seed',        42);
end

fprintf('Puntos: %d densidades T = [%s], esquema unico "%s"\n', ...
    numel(pts), num2str(Tlist), scheme.name);

%% 3. Precalentar la cache P.618 en el CLIENTE (fuera de los cronometros)
%  atm_loss_dB construye su tabla ITU-R P.618 una vez por proceso (~90 s). Si no
%  se precalentara, el modo secuencial pagaria esa construccion y el "speedup"
%  medido seria un artefacto. Los workers se precalientan dentro de
%  run_sweep_points (opt.warmup), tambien fuera del cronometro.
fprintf('\n[warmup] construyendo la tabla P.618 en el cliente...\n');
tw = tic;  atm_loss_dB(cfg, 45);  fprintf('[warmup] listo (%.1f s)\n', toc(tw));

%% 4. (a) BARRIDO SECUENCIAL
[Rseq, Iseq] = run_sweep_points(cfg, pts, 0, struct('warmup',false));

%% 5. (b) BARRIDO PARALELO
[Rpar, Ipar] = run_sweep_points(cfg, pts, [], struct('warmup',true));

%% 6. COMPARACION
fprintf('\n===================== COMPARACION =====================\n');

% OUT.runinfo guarda el tiempo de calculo del punto: es NO determinista por
% definicion y se excluye. Todo lo demas (KPIs, viabilidad, diagnostico,
% definicion de casos, tamanos) entra en la comparacion.
Aseq = cellfun(@(r) rmfield(r,'runinfo'), Rseq, 'UniformOutput', false);
Apar = cellfun(@(r) rmfield(r,'runinfo'), Rpar, 'UniformOutput', false);

eq = isequaln(Aseq, Apar);
fprintf('isequaln(secuencial, paralelo)  : %s\n', ternary(eq,'SI','NO'));

[md, where] = deep_maxdiff(Aseq, Apar, 'results');
fprintf('max|diferencia| numerica        : %.3g   (campo: %s)\n', md, where);

tol = 1e-12;
if md > tol || ~isfinite(md)
    fprintf('\n*** FALLO DE INVARIANCIA ***\n');
    fprintf('El campo %s difiere en %.6g (> %.0e).\n', where, md, tol);
    error('test_parfor_invariance:mismatch', ...
        ['El barrido en paralelo NO reproduce el secuencial (campo %s, dif %.6g). ' ...
         'PARAR: revisar el estado compartido antes de lanzar E3.'], where, md);
end

%% 7. Tabla de resultados (para ver que el barrido calcula algo con sentido)
fprintf('\n%-10s %6s %8s %10s %11s %11s %12s\n', ...
    'punto','N','cob.','SINR_p5_B','R_p5_B[Mb]','margen[dB]','veredicto');
fprintf('%s\n', repmat('-',1,74));
for i = 1:numel(Rseq)
    r = Rseq{i};  V = r.K{1}.viab;
    fprintf('%-10s %6d %8.2f %10.2f %11.1f %11.2f %12s\n', ...
        r.label, r.N, r.diag.covFrac, V.SINR_edge_p5, r.K{1}.edge.R_p5_Mbps, ...
        V.margin_floor_dB, upper(V.verdict));
end
fprintf('%s\n', repmat('-',1,74));

%% 8. Veredicto y rendimiento
sp = Iseq.elapsed_s / max(Ipar.elapsed_s, eps);
fprintf('\n===================== VEREDICTO =====================\n');
fprintf('INVARIANCIA: OK  (isequaln = %s, max|dif| = %.3g <= %.0e)\n', ...
    ternary(eq,'SI','NO'), md, tol);
fprintf('  secuencial : %7.1f s  (%s)\n', Iseq.elapsed_s, Iseq.mode);
fprintf('  paralelo   : %7.1f s  (%s)\n', Ipar.elapsed_s, Ipar.mode);
fprintf('  SPEEDUP    : x%.2f con %d worker(s)\n', sp, Ipar.nWorkers_eff);
if Ipar.nWorkers_eff == 0
    fprintf(['  [Tarea 4] No hay pool disponible: parfor ha degradado a for sin error\n' ...
             '            y el resultado es el mismo. Robustez confirmada.\n']);
elseif numel(pts) < Ipar.nWorkers_eff
    fprintf(['  [nota] Con %d puntos y %d workers el speedup esta limitado por el numero\n' ...
             '         de puntos (techo x%d), no por la maquina.\n'], ...
        numel(pts), Ipar.nWorkers_eff, numel(pts));
end
fprintf('  cache P.618: %d tabla(s) distinta(s) en el barrido -> 1 construccion por worker\n', ...
    Ipar.nAtmKeys);
fprintf('=====================================================\n');

% =========================================================================
function [md, where] = deep_maxdiff(a, b, path)
%DEEP_MAXDIFF  Maxima diferencia absoluta recorriendo structs/cells recursivamente.
%   Devuelve Inf (y la ruta) si las estructuras no son comparables: distinto tipo,
%   distinto tamano, distintos campos, o NaN en uno y no en el otro. Asi una
%   discrepancia estructural no puede pasar como "diferencia 0".
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
    both = isnan(x) & isnan(y);        % NaN vs NaN cuenta como igual (isequaln)
    only = xor(isnan(x), isnan(y));    % NaN en uno solo = incomparable
    d(both) = 0;
    d(only) = Inf;
    md = max(d(:));
    return;
end

% Cualquier otro tipo: comparacion estricta
if ~isequaln(a, b), md = Inf;  where = sprintf('%s [tipo %s]', path, class(a)); end
end

% =========================================================================
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
