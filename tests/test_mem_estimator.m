%TEST_MEM_ESTIMATOR  Coherencia del estimador de memoria con lo MEDIDO en FASE A.
%
%  El estimador de memoria decide cuantos workers de parfor caben en la maquina;
%  si se equivoca a la baja, E3 muere por falta de RAM a mitad del barrido. Este
%  test lo contrasta con los picos REALES medidos en run_convergence_studyA con la
%  configuracion convergida (step_km = 4 -> M = 317 usuarios, Nt = 61, 91 haces):
%
%      T = 1584, sin trocear        ->  2.18 GB  (medido)
%      T = 1584, timeBlock = 20     ->  0.71 GB  (medido, ×3.05 de reduccion)
%      T = 4000, sin trocear        ->  5.40 GB  (medido)
%
%  El estimador ANTERIOR (25 B/elem, sin troceado) daba 0.72 GB para la ventana
%  ENTERA a T = 1584, es decir el valor del BLOQUE atribuido a algo 3x mayor. Ese
%  es exactamente el fallo que se corrige.
%
%  NO ejecuta ningun punto del barrido: solo construye la rejilla y el layout, asi
%  que tarda segundos. No toca la fisica ni el troceado.

clear; clc;
fprintf('===== TEST DEL ESTIMADOR DE MEMORIA (coherencia con la FASE A) =====\n\n');

tol_rel = 0.02;                 % 2%: el modelo debe reproducir lo medido, no aproximarlo

%% Configuracion convergida de la FASE A / E3
cfg = config_default();
cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 4;                 % FASE A: convergido
cfg.time.dt          = 120;               % FASE A: Nt = 61 en 2 h
cfg.beams.nRings     = 5;                 % 91 haces (cluster de run_ffr_demo)
cfg.beams.nBeams     = [];

users = build_user_grid(cfg);
Nt    = numel(cfg.time.t0 : cfg.time.dt : cfg.time.t0 + cfg.time.duration);
BL    = build_beam_layout(cfg);
fprintf('Escenario: M = %d usuarios | Nt = %d instantes | nBeams = %d\n\n', ...
    users.M, Nt, BL.nBeams);
if users.M ~= 317 || Nt ~= 61 || BL.nBeams ~= 91
    warning('test_mem_estimator:scenario', ...
        ['El escenario no reproduce el de la FASE A (317/61/91): las referencias ' ...
         'medidas no son comparables. Revisa config_default.']);
end

%% Casos: [T, timeBlock, pico medido en la FASE A]
cases = { ...
    1584, [], 2.18, 'T=1584 sin trocear' ; ...
    1584, 20, 0.71, 'T=1584 timeBlock=20' ; ...
    4000, [], 5.40, 'T=4000 sin trocear' };

fprintf('%-24s %10s %12s %12s %10s\n', 'caso', 'medido[GB]', 'estimado[GB]', 'dif rel', 'veredicto');
ok = true;
for i = 1:size(cases,1)
    cfgC = cfg;
    cfgC.compute.timeBlock = cases{i,2};
    pt   = struct('T', cases{i,1});

    MEM = estimate_sweep_memory(cfgC, {pt}, 1);
    est = MEM.peak_worker_GB;
    ref = cases{i,3};
    rel = (est - ref) / ref;

    pass = abs(rel) <= tol_rel;
    ok   = ok && pass;
    if pass, v = 'OK'; else, v = 'FALLA'; end
    fprintf('%-24s %10.2f %12.2f %11.1f%% %10s\n', cases{i,4}, ref, est, 100*rel, v);
end

%% El fallo concreto que se corrige: el modelo antiguo de 25 B/elem
old25 = users.M*1584*Nt*(3*8+1) / 2^30;                  % geometria sola, ventana entera
newTB = mem_peak_model_GB(users.M, 1584, 20, BL.nBeams); % modelo real, bloque de 20
fprintf('\nContraste con el estimador ANTERIOR (T=1584):\n');
fprintf('  25 B/elem, ventana entera  -> %.2f GB   (era lo que se imprimia; el real es 2.18)\n', old25);
fprintf('  25 B/elem x factor bloque  -> %.2f GB   (el error si solo se anadiera el troceado)\n', old25*20/Nt);
fprintf('  74 B/elem, bloque de 20    -> %.2f GB   (correcto: coincide con lo medido)\n', newTB);

%% Linealidad en el tamano de bloque (base de la recomendacion de timeBlock)
p20 = mem_peak_model_GB(users.M, 1584, 20, BL.nBeams);
p10 = mem_peak_model_GB(users.M, 1584, 10, BL.nBeams);
linErr = abs(p20 - 2*p10) / p20;
fprintf('\nLinealidad pico(nBlk): |p20 - 2*p10|/p20 = %.3g (debe ser 0)\n', linErr);
ok = ok && linErr < 1e-12;

%% Veredicto
fprintf('\n');
if ok
    fprintf('RESULTADO: el estimador reproduce los picos medidos en la FASE A (<= %.0f%%).\n', 100*tol_rel);
else
    error('test_mem_estimator:mismatch', ...
        'El estimador NO reproduce los picos medidos: revisar mem_peak_model_GB.');
end
