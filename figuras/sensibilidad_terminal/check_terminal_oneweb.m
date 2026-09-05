%% CHECK_TERMINAL_ONEWEB   Sensibilidad de los veredictos al terminal de usuario
%
%  QUE HACE. El perfil de validacion cruzada (config_oneweb) usa a proposito EL
%  MISMO terminal de 60 cm que el perfil principal, aunque el operador real
%  comercialice platos mayores. Este script convierte esa eleccion metodologica en
%  un DATO: reproduce la seccion 3 de run_oneweb_demo (comparativa de los 6
%  esquemas en la mascara nominal) con un terminal de 73 cm y comprueba si algun
%  veredicto de viabilidad cambia.
%
%  NO es un experimento nuevo ni un resultado del trabajo: es una comprobacion de
%  sensibilidad y su unico producto es una tabla.
%
%  QUE NO HACE
%    - no toca ningun fichero del repositorio (config_oneweb, run_oneweb_demo,
%      motor y helpers quedan intactos);
%    - no sobrescribe oneweb_results.mat ni nada de figs_oneweb/;
%    - no relanza el runner completo: solo la seccion 3;
%    - no genera figuras;
%    - la REFERENCIA no se recalcula, se LEE del .mat publicado.
%
%  Salidas (todas en figuras/sensibilidad_terminal/):
%    terminal_check.mat   con cfg, cfg2 y OUT2
%
%  Uso: run('.../figuras/sensibilidad_terminal/check_terminal_oneweb.m')

HERE     = fileparts(mfilename('fullpath'));
PROJROOT = fileparts(fileparts(HERE));
addpath(PROJROOT);

fprintf('\n=======================================================================\n');
fprintf('  SENSIBILIDAD AL TERMINAL DE USUARIO -- escenario OneWeb\n');
fprintf('=======================================================================\n');

%% ============ PASO 1. Referencia publicada (NO se recalcula) =============
REF = load(fullfile(PROJROOT,'oneweb_results.mat'), ...
           'OUT','CASES','nm','cfg','ALPHA_STAR','TAU_Q');
CASES = REF.CASES;      % NO se redefinen: los 6 esquemas y sus parametros son
nm    = REF.nm;         % exactamente los publicados
nC    = numel(CASES);
OUT1  = REF.OUT;

fprintf('\n[PASO 1] Referencia leida de oneweb_results.mat\n');
fprintf('  esquemas   : %s\n', strjoin(nm,' | '));
fprintf('  ALPHA_STAR = %.2f | TAU_Q = %g  (se toman del .mat, no se ajustan)\n', ...
        REF.ALPHA_STAR, REF.TAU_Q);

%% ============ PASO 2. Reconstruir la configuracion =======================
%  Literal de run_oneweb_demo.m seccion 1 (lineas 38-44).
cfg = config_oneweb('hex');
cfg.time.dt           = 120;
cfg.time.duration     = 7200;
cfg.ffr.tau_mode      = 'quantile';
cfg.ffr.tau_q         = REF.TAU_Q;
cfg.ffr.sched         = 'share';
cfg.compute.timeBlock = 10;

%% ============ PASO 3. Coherencia con el cfg publicado (ABORTA) ===========
%  Si difieren, el .mat vendria de otra configuracion y la comparacion no valdria.
fprintf('\n[PASO 3] Coherencia de cfg reconstruida vs cfg guardada en el .mat\n');
campos = {'constellations','geom','radio','beams','ground','time','ffr'};
difs = {};
for k = 1:numel(campos)
    f = campos{k};
    if ~isfield(cfg,f) || ~isfield(REF.cfg,f)
        difs{end+1} = sprintf('%s (ausente en una de las dos)', f); %#ok<SAGROW>
        continue
    end
    A = cfg.(f);  B = REF.cfg.(f);
    if isequaln(A,B)
        fprintf('  cfg.%-15s IDENTICO\n', f);
    else
        % detallar hoja a hoja para poder reportarlo
        if isstruct(A) && isstruct(B) && isscalar(A) && isscalar(B)
            fn = union(fieldnames(A), fieldnames(B));
            sub = {};
            for j = 1:numel(fn)
                hasA = isfield(A,fn{j});  hasB = isfield(B,fn{j});
                if ~hasA
                    sub{end+1} = sprintf('%s (solo en la guardada)', fn{j}); %#ok<SAGROW>
                elseif ~hasB
                    sub{end+1} = sprintf('%s (solo en la reconstruida)', fn{j}); %#ok<SAGROW>
                elseif ~isequaln(A.(fn{j}), B.(fn{j}))
                    sub{end+1} = fn{j}; %#ok<SAGROW>
                end
            end
            fprintf('  cfg.%-15s DIFIERE en: %s\n', f, strjoin(sub,', '));
            difs{end+1} = sprintf('%s -> %s', f, strjoin(sub,', ')); %#ok<SAGROW>
        else
            fprintf('  cfg.%-15s DIFIERE\n', f);
            difs{end+1} = f; %#ok<SAGROW>
        end
    end
end
if ~isempty(difs)
    error('check_terminal_oneweb:cfgMismatch', ...
        ['La cfg reconstruida NO coincide con la guardada en oneweb_results.mat.\n' ...
         'Campos afectados: %s\nLa comparacion no seria valida: PARAR.'], ...
        strjoin(difs, ' | '));
end
fprintf('  -> coherencia OK en los %d bloques comprobados.\n', numel(campos));

%% ============ PASO 4. Variante con el terminal real (73 cm) ==============
lambda = 0.299792458 / cfg.radio.freq_GHz;      % m
D_ref  = 0.60;                                  % m, terminal del perfil principal
D_m    = 0.73;                                  % m, terminal real del operador
DL     = D_m / lambda;
DL_ref = D_ref / lambda;
Tsys   = cfg.radio.Tsys_K;                      % se CONSERVA: cambiar la apertura no
                                                % altera la cadena de ruido del receptor

cfg2 = cfg;
cfg2.radio.Grx_max_dBi = 20*log10(DL) + 7.7;
cfg2.radio.GT_dBK      = cfg2.radio.Grx_max_dBi - 10*log10(Tsys);
cfg2.radio.Tsys_K      = 10^((cfg2.radio.Grx_max_dBi - cfg2.radio.GT_dBK)/10);
assert_cfg_coherent(cfg2);

dG = cfg2.radio.Grx_max_dBi - cfg.radio.Grx_max_dBi;

fprintf('\n[PASO 4] Terminales comparados\n');
fprintf('  CRITERIO: la figura de merito G/T es el parametro MAESTRO en este perfil,\n');
fprintf('  igual que en el principal. Se mantiene Tsys porque un cambio de apertura no\n');
fprintf('  altera la cadena de ruido del receptor, y el incremento de ganancia se\n');
fprintf('  traslada INTEGRO a la figura de merito. term_efficiency NO se toca.\n\n');
fprintf('  %-26s %14s %14s %10s\n','magnitud','60 cm (ref)','73 cm (real)','dif');
fprintf('  %-26s %14.4f %14.4f %10s\n','diametro [m]', D_ref, D_m, '');
fprintf('  %-26s %14.4f %14.4f %10s\n','D/lambda geometrico', DL_ref, DL, '');
fprintf('  %-26s %14.4f %14.4f %+10.4f\n','Grx_max [dBi]', ...
        cfg.radio.Grx_max_dBi, cfg2.radio.Grx_max_dBi, dG);
fprintf('  %-26s %14.4f %14.4f %+10.4f\n','G/T [dB/K]', ...
        cfg.radio.GT_dBK, cfg2.radio.GT_dBK, cfg2.radio.GT_dBK - cfg.radio.GT_dBK);
fprintf('  %-26s %14.4f %14.4f %+10.4f\n','Tsys [K]', ...
        cfg.radio.Tsys_K, cfg2.radio.Tsys_K, cfg2.radio.Tsys_K - cfg.radio.Tsys_K);
fprintf('  %-26s %14.4f %14.4f\n','term_efficiency [-]', ...
        cfg.radio.term_efficiency, cfg2.radio.term_efficiency);
fprintf('  lambda = %.6f m\n', lambda);

% Coherencia del metodo: la MISMA formula aplicada al plato de 60 cm debe
% reproducir el Grx_max que trae config_oneweb. Si no, el criterio del paso 4 no
% seria "el mismo metodo" y habria que reportarlo.
G60_formula = 20*log10(DL_ref) + 7.7;
fprintf('\n  [coherencia del metodo] 20log10(D/lambda)+7.7 con D=0.60 m da %.4f dBi;\n', G60_formula);
fprintf('  config_oneweb declara %.4f dBi -> desviacion %.4f dB.\n', ...
        cfg.radio.Grx_max_dBi, G60_formula - cfg.radio.Grx_max_dBi);
if abs(G60_formula - cfg.radio.Grx_max_dBi) > 0.05
    warning('check_terminal_oneweb:metodo', ...
        ['La formula no reproduce el Grx_max de 60 cm (desviacion %.3f dB): ' ...
         'el terminal nuevo NO estaria derivado con el mismo metodo.'], ...
        G60_formula - cfg.radio.Grx_max_dBi);
end

%% ---- Cortes de la envolvente S.1428-1 y aviso de regimen ----------------
%  term_gain_dB deriva D/lambda de (Grx_max, term_efficiency), no del diametro:
%  se reportan los dos para que se vea si coinciden.
env = @(c) struct( ...
    'DL',    sqrt(10^(c.radio.Grx_max_dBi/10) / (c.radio.term_efficiency*pi^2)), ...
    'G1',    29 - 25*log10(95/sqrt(10^(c.radio.Grx_max_dBi/10)/(c.radio.term_efficiency*pi^2))), ...
    'Gmax',  c.radio.Grx_max_dBi);
e1 = env(cfg);   e2 = env(cfg2);
e1.phi_m = (20/e1.DL)*sqrt(max(e1.Gmax - e1.G1,0));  e1.phi_r = 95/e1.DL;
e2.phi_m = (20/e2.DL)*sqrt(max(e2.Gmax - e2.G1,0));  e2.phi_r = 95/e2.DL;

fprintf('\n  --- cortes de la envolvente ITU-R S.1428-1 ---\n');
fprintf('  %-30s %14s %14s\n','', '60 cm (ref)','73 cm (real)');
fprintf('  %-30s %14.4f %14.4f\n','D/lambda derivado en term_gain_dB', e1.DL, e2.DL);
fprintf('  %-30s %14.4f %14.4f\n','G1, meseta [dBi]', e1.G1, e2.G1);
fprintf('  %-30s %14.4f %14.4f\n','phi_m, fin lob. principal [deg]', e1.phi_m, e2.phi_m);
fprintf('  %-30s %14.4f %14.4f\n','phi_r, inicio 29-25log [deg]', e1.phi_r, e2.phi_r);
fprintf('  (los escalones de cola, 33.1 y 80 deg, son fijos en la Recomendacion)\n');

phis = [0 1 2 5 10 33.1 45 80 120 170];
fprintf('\n  ganancia [dBi] en angulos de muestra:\n');
fprintf('  %8s', 'phi[deg]'); fprintf(' %8.1f', phis); fprintf('\n');
fprintf('  %8s', '60 cm');    fprintf(' %8.2f', term_gain_dB(cfg,  phis)); fprintf('\n');
fprintf('  %8s', '73 cm');    fprintf(' %8.2f', term_gain_dB(cfg2, phis)); fprintf('\n');

% --- Aviso de regimen: debe DESAPARECER al pasar de D/L~23.9 a ~29.2 -----
fprintf('\n  --- aviso de regimen de term_gain_dB ---\n');
clear term_gain_dB; lastwarn('');
term_gain_dB(cfg, 0);
[w1, id1] = lastwarn;
clear term_gain_dB; lastwarn('');
term_gain_dB(cfg2, 0);
[w2, id2] = lastwarn;
fprintf('  60 cm: %s\n', tern_str(isempty(w1), 'SIN aviso', sprintf('AVISO [%s]', id1)));
if ~isempty(w1), fprintf('         %s\n', w1); end
fprintf('  73 cm: %s\n', tern_str(isempty(w2), 'SIN aviso', sprintf('AVISO [%s]', id2)));
if ~isempty(w2), fprintf('         %s\n', w2); end
warn60 = w1;  warn73 = w2;

%% ============ PASO 5. Ejecutar la seccion 3 con el terminal real =========
fprintf('\n[PASO 5] Ejecutando la comparativa de esquemas con el terminal de 73 cm\n');
fprintf('  (misma llamada que run_oneweb_demo seccion 3; CASES y nm salen del .mat)\n\n');
pt = struct('cases',{CASES}, 'elMin_valid',70, 'label','OneWeb terminal 73 cm', ...
            'verbose',true, 'idx',1);
tExec = tic;
OUT2  = run_one_density(cfg2, pt);
tRun  = toc(tExec);
fprintf('\n  [tiempo] %.1f s (%.2f min)\n', tRun, tRun/60);

%% ============ PASO 6. Tabla comparativa =================================
fprintf('\n=======================================================================\n');
fprintf('  TABLA COMPARATIVA -- terminal 60 cm (publicado) vs 73 cm\n');
fprintf('=======================================================================\n');
fprintf('%-14s | %17s | %15s | %13s | %19s\n', ...
        'esquema','SINR_edge_p5 [dB]','R_p5 borde[Mbps]','cobertura','veredicto');
fprintf('%-14s | %8s %8s | %7s %7s | %6s %6s | %9s %9s\n', ...
        '','60cm','73cm','60cm','73cm','60cm','73cm','60cm','73cm');
fprintf('%s\n', repmat('-',1,92));

dSINR = nan(1,nC);  changed = false(1,nC);
for c = 1:nC
    a = OUT1.K{c};  b = OUT2.K{c};
    dSINR(c)  = b.viab.SINR_edge_p5 - a.viab.SINR_edge_p5;
    changed(c)= ~strcmp(a.viab.verdict, b.viab.verdict);
    fprintf('%-14s | %+8.3f %+8.3f | %7.2f %7.2f | %6.2f %6.2f | %9s %9s %s\n', ...
        nm{c}, a.viab.SINR_edge_p5, b.viab.SINR_edge_p5, ...
        a.edge.R_p5_Mbps, b.edge.R_p5_Mbps, ...
        a.viab.coverage, b.viab.coverage, ...
        a.viab.verdict, b.viab.verdict, tern_str(changed(c),'  <== CAMBIA',''));
end
fprintf('%s\n', repmat('-',1,92));
fprintf('%-14s | %+8s %+8.3f |\n','diferencia SINR','', max(abs(dSINR)));

%  MARGENES del caso de REFERENCIA. El veredicto de 3 niveles tiene DOS fronteras
%  (suelo QPSK y umbral de servicio), asi que un desplazamiento de SINR puede mover
%  el veredicto cruzando CUALQUIERA de las dos: se comprueban las dos, no solo el
%  umbral, porque con solo margin_th_dB el argumento quedaria incompleto.
fprintf('\n  margenes del caso de REFERENCIA (las DOS fronteras del veredicto):\n');
fprintf('  gamma_floor = %+.2f dB | gamma_th = %+.2f dB\n', ...
        OUT1.K{1}.viab.gamma_floor_dB, OUT1.K{1}.viab.gamma_th_dB);
fprintf('  %-14s %14s %14s %12s   %s\n', ...
        'esquema','margen_th[dB]','margen_suelo[dB]','|dSINR|[dB]','cabe en la frontera mas cercana?');
for c = 1:nC
    mt = OUT1.K{c}.viab.margin_th_dB;
    mf = OUT1.K{c}.viab.margin_floor_dB;
    dmin = min(abs(mt), abs(mf));                 % frontera MAS CERCANA
    fprintf('  %-14s %+14.3f %+14.3f %12.3f   %s\n', nm{c}, mt, mf, abs(dSINR(c)), ...
        tern_str(abs(dSINR(c)) < dmin, ...
            sprintf('SI (holgura x%.0f)', dmin/max(abs(dSINR(c)),eps)), ...
            'NO -- podria mover el veredicto'));
end

fprintf('\n  VEREDICTOS QUE CAMBIAN: %d de %d', nnz(changed), nC);
if any(changed)
    fprintf('  -> %s\n', strjoin(nm(changed),', '));
    fprintf('  *** HALLAZGO: hay que decidirlo a mano. NO se re-ejecuta. ***\n');
else
    fprintf('  (ninguno)\n');
end
fprintf('  Mayor diferencia de SINR observada: %.4f dB (%s)\n', ...
        max(abs(dSINR)), nm{find(abs(dSINR)==max(abs(dSINR)),1)});

%% ============ PASO 7. Desglose de potencias (lo que este accesible) ======
fprintf('\n--- desglose de potencias (reuse1 y ffr(3,a*)) ---\n');
fprintf(['NOTA: run_one_density acumula SOLO SINR_dB, CN_dB y penalty_dB\n' ...
         '(lineas 266-268) y compute_kpis guarda penalty_mean/p5/p95, pero NO\n' ...
         'propaga I_intra ni I_inter por separado. El desglose intra/inter NO es\n' ...
         'accesible en la salida y se OMITE en lugar de recalcularlo.\n' ...
         'Lo que si se puede reportar sin recalcular nada:\n' ...
         '  penalty = C/N - SINR = 10log10(1 + I_total/N)  -> I_total/N\n' ...
         '  C/N media = SINR_media + penalty_media (misma muestra, luego exacto)\n']);
sel = [1 4];
fprintf('\n  %-12s %-6s %10s %10s %10s %12s\n', ...
        'esquema','term','SINR[dB]','penal[dB]','C/N[dB]','I_tot/N[dB]');
for c = sel
    for t = 1:2
        if t==1, K = OUT1.K{c}; lab='60cm'; else, K = OUT2.K{c}; lab='73cm'; end
        sinrm = K.all.SINR_mean;  penm = K.penalty_mean;
        cn    = sinrm + penm;
        inr   = 10*log10(max(10^(penm/10) - 1, realmin));
        fprintf('  %-12s %-6s %10.3f %10.3f %10.3f %12.3f\n', nm{c}, lab, sinrm, penm, cn, inr);
    end
end
fprintf(['\n  LECTURA: si la ganancia del terminal se cancela en C/I_intra y solo\n' ...
         '  sobrevive frente al ruido termico, C/N debe subir ~%.2f dB (el dG del\n' ...
         '  terminal) y la SINR quedarse practicamente igual.\n'], dG);

%% ============ PASO 8. Guardar ============================================
save(fullfile(HERE,'terminal_check.mat'), 'cfg','cfg2','OUT2');
fprintf('\nGuardado %s\n', fullfile(HERE,'terminal_check.mat'));
fprintf('(oneweb_results.mat y figs_oneweb/ NO se han tocado)\n');

%% ------------------------------------------------------------------------
function s = tern_str(c, a, b)
if c, s = a; else, s = b; end
end
