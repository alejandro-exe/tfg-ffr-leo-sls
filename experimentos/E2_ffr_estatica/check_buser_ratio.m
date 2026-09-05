%% CHECK_BUSER_RATIO  Verificacion de la derivacion de B_user en reuseD vs FFR.
%  SCRIPT AUXILIAR de comprobacion. NO forma parte del corpus de experimentos:
%  no guarda .mat, no genera figuras y no toca el motor.
%
%  PREGUNTA. En la tabla de E2 reuseD y ffr dan la MISMA
%  SINR_edge_p5 pero distinto R_p5 de borde. Como R = B_user*log2(1+SINR), toda la
%  diferencia debe estar en B_user. La derivacion a verificar es:
%
%      reuseD (una sola clase, f_clase = 1):   B_user = B / (Delta*n0)
%      ffr, usuario de BORDE:                  B_user = (1-alpha)*B / (Delta*n0*f_borde)
%      cociente predicho = (1-alpha)/f_borde = 0.60/0.50 = 1.20
%
%  El punto flojo es si en reuseD el codigo divide tambien por la fraccion de
%  clase (entonces el cociente seria 0.60). Esto lo mide.
%
%  CONFIGURACION BARATA A PROPOSITO: el cociente es aritmetica de reparto de banda
%  (Sec.7 de ffr_allocate) y no depende de la constelacion ni de la densidad, asi
%  que NO se reproduce E2 a T=1584. T=66, paso 4 km, Nt=31.
%
%  CLAVE METODOLOGICA: los cuatro casos se evaluan sobre el MISMO ctx y con la
%  MISMA SINR de referencia (precalculada una vez), para que n0, tau y la
%  clasificacion centro/borde sean identicos entre esquemas.

clear; clc;
tRun = tic;

%% 1. Configuracion barata (construida aqui; config_default NO se toca)
cfg = config_default();

cfg.constellations(1).T = 66;          % constelacion reducida: el cociente no depende de T
cfg.constellations(1).P = 6;

cfg.ground.mode      = 'localgrid';
cfg.ground.radius_km = 40;
cfg.ground.step_km   = 4;              % perfil de FASE A

cfg.time.dt       = 240;               % s
cfg.time.duration = 7200;              % s (2 h) -> Nt = 31

cfg.beams.nRings = 5;                  % 91 haces

cfg.ffr.classifier = 'sinr';
cfg.ffr.tau_mode   = 'quantile';
cfg.ffr.tau_q      = 50;               % 50% de CENTRO -> f_borde = 0.50
cfg.ffr.sched      = 'share';
cfg.ffr.adaptive   = false;            % alpha constante (FFR estatica)

cfg.compute.timeBlock = [];            % sin trocear: el punto es diminuto

fprintf('===== CHECK B_user: reuseD vs FFR (script desechable) =====\n');

% Warmup de la tabla P.618 fuera del cronometro del experimento (coste fijo de la
% configuracion atmosferica, cacheado en un persistent).
tW = tic;  atm_loss_dB(cfg, 45);  tWarm = toc(tW);
fprintf('[warmup] tabla P.618: %.1f s\n', tWarm);
tExp = tic;

%% 2. Contexto UNICO compartido por los cuatro casos
[ctx, sats, users, BL] = build_ctx_blocked(cfg, true);
fprintf('M=%d usuarios | Nt=%d | N=%d satelites | %d haces | s=%.4f deg\n', ...
    users.M, ctx.Nt, sats.N, BL.nBeams, BL.spacing_deg);

% SINR de referencia (reuso-1) precalculada UNA vez -> misma tau y misma
% clasificacion centro/borde en los cuatro casos (igual que en run_ffr_demo).
cfgR = cfg;  cfgR.ffr.scheme = 'reuse1';  cfgR.ffr.Delta = 1;
AREF = ffr_allocate(cfgR, ctx, false);
ctx.SINR_ref_dB = AREF.SINR_ref_dB;

%% 3. Los cuatro casos, sobre el MISMO ctx
cases = { ...
    struct('scheme','reuseD','Delta',3,'alpha',NaN), ...
    struct('scheme','ffr',   'Delta',3,'alpha',0.40), ...
    struct('scheme','reuseD','Delta',4,'alpha',NaN), ...
    struct('scheme','ffr',   'Delta',4,'alpha',0.40) };
nC = numel(cases);
AA = cell(1,nC);
fB = nan(1,nC);  fBall = nan(1,nC);
mE = nan(1,nC);  sE = nan(1,nC);
mC = nan(1,nC);  sC = nan(1,nC);

for c = 1:nC
    cc = cases{c};
    cfg.ffr.scheme = cc.scheme;
    cfg.ffr.Delta  = cc.Delta;
    if ~isnan(cc.alpha), cfg.ffr.alpha = cc.alpha; end

    A = ffr_allocate(cfg, ctx, false);
    AA{c} = A;

    isC = A.isCenter;  isE = ctx.cov & ~isC;      % borde = cubierto y no centro
    Bu  = A.B_user_Hz;

    fBall(c) = mean(~isC(:));                     % definicion literal (todas las muestras)
    fB(c)    = sum(isE(:)) / max(sum(ctx.cov(:)),1);  % la que usa el codigo (solo cubiertos)
    mE(c) = mean(Bu(isE),'omitnan')/1e6;  sE(c) = std(Bu(isE),'omitnan')/1e6;
    mC(c) = mean(Bu(isC),'omitnan')/1e6;  sC(c) = std(Bu(isC),'omitnan')/1e6;
end

% Comprobacion: la clasificacion debe ser IDENTICA en los cuatro casos
for c = 2:nC
    if ~isequal(AA{c}.isCenter, AA{1}.isCenter)
        error('check_buser_ratio:pop', 'La clasificacion centro/borde difiere entre casos.');
    end
end

%% 4. Salida
fprintf('\n');
fprintf('esquema        Delta  alpha   f_borde   B_user_borde[MHz]   B_user_centro[MHz]\n');
for c = 1:nC
    cc = cases{c};
    if isnan(cc.alpha), as = '  -  '; else, as = sprintf('%.2f', cc.alpha); end
    fprintf('%-14s %3d    %5s   %6.4f      %8.4f            %8.4f\n', ...
        cc.scheme, cc.Delta, as, fB(c), mE(c), mC(c));
end
fprintf('\n');
fprintf('cociente B_user_borde  ffr(3)/reuseD(3) = %.4f\n', mE(2)/mE(1));
fprintf('cociente B_user_borde  ffr(4)/reuseD(4) = %.4f\n', mE(4)/mE(3));

% Identidad B = Bc + Delta*Be sobre los casos FFR (donde Bc/Be se usan de verdad)
B = radio_band_Hz(cfg);
e3 = max(abs(AA{2}.Bc_Hz + AA{2}.Delta*AA{2}.Be_Hz - B));
e4 = max(abs(AA{4}.Bc_Hz + AA{4}.Delta*AA{4}.Be_Hz - B));
fprintf('identidad B = Bc + D*Be   (D=3): %.3g   (D=4): %.3g     <- max|error| en Hz\n', e3, e4);

%% 5. Notas: dispersion dentro de clase y definicion de f_borde
fprintf('\n--- notas (no pedidas en la tabla, pero relevantes) ---\n');
fprintf('desv.tipica de B_user dentro de clase [MHz]:\n');
for c = 1:nC
    fprintf('  %-8s D=%d : borde %.4f (%.2f%% de la media) | centro %.4f (%.2f%%)\n', ...
        cases{c}.scheme, cases{c}.Delta, sE(c), 100*sE(c)/mE(c), sC(c), 100*sC(c)/mC(c));
end
fprintf('n0 (usuarios/celda): media %.4f, rango [%.4f %.4f] sobre instantes con cobertura\n', ...
    mean(AA{1}.load_n0(sum(ctx.cov,1)>0)), min(AA{1}.load_n0(sum(ctx.cov,1)>0)), ...
    max(AA{1}.load_n0(sum(ctx.cov,1)>0)));
fprintf('f_borde sobre muestras CUBIERTAS = %.4f (es la que usa ffr_allocate)\n', fB(1));
fprintf('f_borde con la definicion literal mean(~isCenter) sobre TODAS las muestras = %.4f\n', fBall(1));
fprintf('cobertura = %.4f (%d de %d muestras)\n', ...
    mean(ctx.cov(:)), sum(ctx.cov(:)), numel(ctx.cov));

% Prediccion analitica frente a lo medido
fprintf('\nprediccion (1-alpha)/f_borde = %.4f/%.4f = %.4f\n', ...
    1-0.40, fB(1), (1-0.40)/fB(1));

fprintf('\n[tiempo] experimento %.1f s (+ %.1f s de tabla P.618) | total %.1f s\n', ...
    toc(tExp), tWarm, toc(tRun));
