function export_h2_criterion()
%EXPORT_H2_CRITERION  Consolida el criterio de rentabilidad del interior de la FFR.
%
%   POST-PROCESO PURO: NO simula, NO toca el motor, NO genera figuras y NO
%   escribe fuera de h2_rerun/. Lee tres .mat ya publicados y los reduce a UNA
%   tabla plana con el criterio de H2 y sus DOS factores por separado.
%
%   EL CRITERIO (experimento E5). Con sched='share' el agregado por
%   celda es  R = Bc*SE_c + Be*SE_b,  y sustituyendo la particion estricta
%   Bc = alpha*B, Be = (1-alpha)*B/Delta y derivando respecto de alpha:
%
%        dR/dalpha = B*(SE_c - SE_b/Delta)   =>   el INTERIOR renta  <=>  crit > 1,
%        con  crit = Delta * SE_centro / SE_borde.
%
%   Es el criterio de validez de H2: si crit < 1 en TODA la ventana de operacion,
%   no hay ninguna geometria en la que subir alpha compense, y por tanto la FFR
%   adaptativa no puede ganar POR CONSTRUCCION (resultado ESTRUCTURAL, no un
%   problema de ajuste de anclas).
%
%   POR QUE SE EMITEN LOS DOS FACTORES POR SEPARADO. El cociente solo dice si el
%   interior renta; SE_c y SE_b dicen POR QUE. Un crit < 1 puede venir de un
%   centro pobre (el interior va en reuso-1, sin proteccion) o de un borde muy
%   bueno (reuso Delta), y en la defensa hay que poder distinguirlo. Por eso la
%   tabla lleva SE_c, SE_b y crit, no solo crit.
%
%   FUENTES (los tres .mat se leen tal cual, no se regeneran):
%     e5_adaptive_results.mat  ratio_bin/ratio_all: Starlink T=1584, estratificado
%                              por tramo de elevacion + global (E5, cap. 7).
%     e3_density_sweep.mat     results{i}.diag.crit_ratio_ffr: eje de densidad
%                              orbital T = 66..4000 (v7.3).
%     oneweb_results.mat       OUT/OUTR.diag.crit_ratio_ffr: validacion cruzada de
%                              arquitectura (hex 19 haces / rect 4x4 16 haces).
%                              RE (barrido de mascara de elevacion, 6 puntos) y
%                              RD (barrido de anillos con la HUELLA FIJA, 4 puntos).
%
%   POR QUE RE Y RD (anadidos en esta revision). Antes solo se leian OUT y OUTR, y
%   en RE/RD vive la MITAD de los cruces del criterio de todo el corpus: sin ellos
%   el apartado "LO QUE CRUZA 1" informaba de menos casos de los que hay. RD(1) se
%   SALTA porque nRings=2 ES la configuracion nominal y duplica exactamente la fila
%   de OUT (se comprueba con un assert; ver Sec. 4c).
%
%   SALIDA:  h2_rerun/results/h2_criterion.mat  ->  tabla H2 + CHK + META.
%   LOG   :  h2_rerun/logs/export_h2_criterion.txt (ademas de la consola).
%
%   VERIFICACION BLOQUEANTE: el script compara las cifras de control contra los
%   valores publicados y ABORTA SIN GUARDAR si alguna no coincide. Es el ancla que
%   impide que un .mat regenerado a espaldas de este post-proceso pase por bueno.
%   DOS TOLERANCIAS, y el log declara cual usa cada comprobacion:
%     TOL9 = 1e-9  para el criterio (crit) de las 24 configuraciones cuyo valor de
%                  referencia esta transcrito a PRECISION COMPLETA;
%     TOL  = 1e-4  para las referencias heredadas a 6 decimales (E5) y para los
%                  controles estructurales (d) y (e), cuyo valor medido es ~1e-16.
%   Una referencia escrita a 6 decimales NO puede comprobarse a 1e-9: su propio
%   redondeo vale ya 5e-7. Por eso subir la tolerancia exige subir la precision de
%   la referencia, y es lo que se ha hecho con las cifras de crit.
%
%   Uso:  export_h2_criterion      (desde la raiz del proyecto o desde h2_rerun/)
%
%   TFG UC3M - Viabilidad de FFR adaptativa en constelaciones LEO.

%% 0. Rutas y log -----------------------------------------------------------
% Las rutas se derivan de la ubicacion DEL PROPIO FICHERO, no del directorio de
% trabajo: asi el script da lo mismo lanzado desde la raiz del proyecto o desde
% h2_rerun/, y no depende de que el path de MATLAB este en un estado concreto.
H2DIR   = fileparts(mfilename('fullpath'));
ROOTDIR = fileparts(H2DIR);                       % raiz del proyecto
RESDIR  = fullfile(H2DIR, 'results');
LOGDIR  = fullfile(H2DIR, 'logs');
if ~exist(RESDIR,'dir'), mkdir(RESDIR); end
if ~exist(LOGDIR,'dir'), mkdir(LOGDIR); end

LOGF = fullfile(LOGDIR, 'export_h2_criterion.txt');
diary off;
if exist(LOGF,'file'), delete(LOGF); end          % log limpio en cada ejecucion
diary(LOGF);
% El log se cierra pase lo que pase (incluida una salida por error()): al ser una
% FUNCION, el objeto de limpieza se destruye al deshacerse su workspace.
cleanupDiary = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('=======================================================================\n');
fprintf('  CRITERIO DE RENTABILIDAD DEL INTERIOR:  Delta*SE_centro/SE_borde > 1\n');
fprintf('  Consolidacion de las TRES fuentes publicadas (post-proceso, sin simular)\n');
fprintf('=======================================================================\n');
stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
fprintf('fecha        : %s\n', stamp);
fprintf('proyecto     : %s\n', ROOTDIR);
fprintf('salida       : %s\n', fullfile(RESDIR,'h2_criterion.mat'));
fprintf('log          : %s\n', LOGF);

TOL  = 1e-4;                                      % <-- DECISION (referencias a 6 decimales)
TOL9 = 1e-9;                                      % <-- DECISION (referencias a precision completa)
fprintf('tolerancias  : crit a precision completa %.0e | referencias heredadas %.0e\n', TOL9, TOL);

%% 1. Comprobacion de que las tres fuentes estan ---------------------------
SRC = {'e5_adaptive_results.mat', 'e3_density_sweep.mat', 'oneweb_results.mat'};
fprintf('\n--- FUENTES ---\n');
for i = 1:numel(SRC)
    f = fullfile(ROOTDIR, SRC{i});
    if ~exist(f,'file')
        error('export_h2_criterion:falta', ...
            'No se encuentra la fuente %s en %s. Este script NO la regenera.', SRC{i}, ROOTDIR);
    end
    d = dir(f);
    fprintf('  %-26s %8.2f MB   %s\n', SRC{i}, d.bytes/2^20, ...
        datestr(d.datenum,'yyyy-mm-dd HH:MM'));
end

% Acumuladores de la tabla plana (una fila = una configuracion evaluada).
%   'variante' = valor del eje barrido (mascara en grados o numero de anillos), NaN
%                en las filas que no pertenecen a un barrido de OneWeb.
%   'eje'      = etiqueta del PUNTO del barrido SIN la politica. Es la parte de
%                'caso' que identifica la configuracion fisica: se usa para la clave
%                de deduplicacion (Sec. 6), donde la estatica y la adaptativa del
%                mismo Delta DEBEN colapsar y dos puntos distintos del mismo eje NO.
%                No se emite a la tabla; su informacion ya esta en caso+variante.
ROW = struct('fuente',{{}}, 'arquitectura',{{}}, 'caso',{{}}, 'variante',[], ...
             'Delta',[], 'T',[], 'nBeams',[], 'elev_mean',[], 'n_muestras',[], ...
             'SE_c',[], 'SE_b',[], 'crit',[], 'eje',{{}});

CHK  = struct();                                  % desviaciones medidas
fail = {};                                        % comprobaciones fuera de tolerancia

%% ========================================================================
%  2. FUENTE 1 - E5: Starlink T=1584, estratificado por elevacion + global
%  ========================================================================
fprintf('\n=======================================================================\n');
fprintf('  1. E5 | Starlink T=1584/72 | estratificado por tramo de elevacion\n');
fprintf('=======================================================================\n');

E5 = load(fullfile(ROOTDIR,'e5_adaptive_results.mat'), ...
          'ratio_bin','ratio_all','DELTA','bin','binName','cases','KK','BL','cfg','ctxLite');

nB_E5 = E5.BL.nBeams;
T_E5  = E5.cfg.constellations(1).T;
D_E5  = E5.DELTA;
iST   = 2;                                        % indice de 'ffr-estatica' en E5
if ~strcmp(E5.cases{iST}.name, 'ffr-estatica')
    error('export_h2_criterion:iST', ...
        'El caso %d de E5 es "%s" y se esperaba "ffr-estatica": el .mat cambio de orden.', ...
        iST, E5.cases{iST}.name);
end

% Tramos de elevacion: EXACTAMENTE los de run_e5_adaptive (Sec. 6). Se replican
% aqui solo para recomputar n y la elevacion media desde ctxLite. El COCIENTE no
% se recalcula: se toma ratio_bin del .mat, que es el numero publicado.
edgesEl = [25 35 60 90];
elevG   = E5.ctxLite.elev_deg;
covG    = E5.ctxLite.cov;

fprintf('%-14s %10s %12s %12s %12s %14s %9s\n', ...
    'tramo','muestras','elev media','SE_centro','SE_borde','D*SEc/SEb','renta?');
fprintf('%s\n', repmat('-',1,90));

for e = 1:3
    mk  = covG & elevG >= edgesEl(e) & elevG < edgesEl(e+1);
    nMk = sum(mk(:));
    if nMk > 0, elMk = mean(elevG(mk)); else, elMk = NaN; end

    % CONTROL: la mascara reconstruida debe reproducir el recuento que el propio
    % E5 guardo en bin.n. Si no coincide, ctxLite y bin son de tandas distintas.
    if nMk ~= E5.bin.n(e)
        error('export_h2_criterion:binN', ...
            ['El tramo %s tiene %d muestras reconstruidas frente a las %d guardadas ' ...
             'en bin.n: ctxLite y bin no son de la misma tanda.'], ...
            E5.binName{e}, nMk, E5.bin.n(e));
    end

    % SE_c y SE_b POR TRAMO NO ESTAN EN EL .mat: run_e5_adaptive los calcula al
    % vuelo (compute_kpis con mascara) para imprimir su tabla y guarda solo el
    % COCIENTE (ratio_bin). Recuperarlos exigiria reejecutar la simulacion, que
    % es justo lo que este post-proceso no hace, asi que van NaN y queda declarado.
    ROW = addrow(ROW, 'E5', 'Starlink-hex', sprintf('bin %s', E5.binName{e}), ...
                 NaN, sprintf('bin %s', E5.binName{e}), ...
                 D_E5, T_E5, nB_E5, elMk, nMk, NaN, NaN, E5.ratio_bin(e));

    fprintf('%-14s %10d %12s %12s %12s %14s %9s\n', E5.binName{e}, nMk, ...
        fmtnum(elMk,'%.2f'), 'n/d', 'n/d', fmtnum(E5.ratio_bin(e),'%.6f'), ...
        rentaTxt(E5.ratio_bin(e)));
end

% Fila GLOBAL: aqui SI hay los dos factores, porque salen de KK{iST}, que si esta
% guardado. Es exactamente el par con el que run_e5_adaptive calculo ratio_all.
SEc_E5 = E5.KK{iST}.center.SE_mean_bpsHz;
SEb_E5 = E5.KK{iST}.edge.SE_mean_bpsHz;
nAll   = sum(covG(:));
elAll  = mean(elevG(covG));
ROW = addrow(ROW, 'E5', 'Starlink-hex', 'global (ventana completa)', NaN, 'global', ...
             D_E5, T_E5, nB_E5, elAll, nAll, SEc_E5, SEb_E5, E5.ratio_all);
fprintf('%-14s %10d %12.2f %12.6f %12.6f %14.6f %9s\n', 'TODO', nAll, elAll, ...
    SEc_E5, SEb_E5, E5.ratio_all, rentaTxt(E5.ratio_all));

% Coherencia interna: ratio_all debe ser Delta*SE_c/SE_b con los SE guardados.
CHK.E5_ratio_all_recompute = abs(D_E5*SEc_E5/SEb_E5 - E5.ratio_all);
fprintf('\n[control] ratio_all recalculado desde KK{%d}: |dif| = %.3e\n', ...
    iST, CHK.E5_ratio_all_recompute);
fprintf(['[nota]    SE_centro y SE_borde POR TRAMO no estan en e5_adaptive_results.mat\n' ...
         '          (run_e5_adaptive guarda solo el cociente ratio_bin). Van como NaN:\n' ...
         '          rellenarlos exigiria reejecutar E5, y este script no simula.\n']);
fprintf(['[nota]    El tramo 25-35 deg esta VACIO (n = 0) y por eso su cociente es NaN:\n' ...
         '          con 1584 satelites el servidor nunca baja de ~49 deg.\n']);

%% ========================================================================
%  3. FUENTE 2 - Eje de densidad orbital (T = 66 .. 4000, P escalado)
%  ========================================================================
fprintf('\n=======================================================================\n');
fprintf('  2. BARRIDO DE DENSIDAD ORBITAL | Starlink, T = 66 .. 4000 (P escalado)\n');
fprintf('=======================================================================\n');

DS   = load(fullfile(ROOTDIR,'e3_density_sweep.mat'), 'results','Tlist','Plist','schemeCases');
nPt  = numel(DS.results);
COLS = [4 5 6];                                   % FFR(D=3) | FFR(D=4) | ADAPT(D=3)

% Guarda de indices: las tres columnas pedidas deben ser esquemas 'ffr'. Si el
% orden de schemeCases cambiara, esto aborta en vez de emitir columnas erroneas.
for k = 1:numel(COLS)
    c = COLS(k);
    if ~strcmpi(DS.schemeCases{c}.scheme, 'ffr')
        error('export_h2_criterion:colDens', ...
            'La columna %d del barrido es "%s" (scheme=%s) y se esperaba un esquema ffr.', ...
            c, DS.schemeCases{c}.name, DS.schemeCases{c}.scheme);
    end
end
nmCols = cell(1,numel(COLS));
for k = 1:numel(COLS), nmCols{k} = DS.schemeCases{COLS(k)}.name; end
fprintf('columnas leidas: %s\n', strjoin(nmCols, ' | '));

fprintf('\n%-6s %-6s %-22s %8s %10s %12s %12s %14s %9s\n', ...
    'T','P','esquema','nBeams','elev med','SE_centro','SE_borde','D*SEc/SEb','renta?');
fprintf('%s\n', repmat('-',1,110));

critDens = nan(nPt, numel(COLS));
for i = 1:nPt
    R = DS.results{i};
    for k = 1:numel(COLS)
        c   = COLS(k);
        SEc = R.K{c}.center.SE_mean_bpsHz;
        SEb = R.K{c}.edge.SE_mean_bpsHz;
        cr  = R.diag.crit_ratio_ffr(c);
        critDens(i,k) = cr;

        % eje = 'nominal': en este barrido el punto lo identifica T, que ya es
        % columna propia de la clave de deduplicacion.
        ROW = addrow(ROW, 'densidad', 'Starlink-hex', R.cases{c}.name, NaN, 'nominal', ...
                     R.cases{c}.Delta, R.T, R.nBeams, R.diag.elev_mean, ...
                     R.K{c}.all.n, SEc, SEb, cr);

        fprintf('%-6d %-6d %-22s %8d %10.2f %12.6f %12.6f %14.6f %9s\n', ...
            R.T, R.P, R.cases{c}.name, R.nBeams, R.diag.elev_mean, SEc, SEb, cr, rentaTxt(cr));
    end
end

%% ========================================================================
%  4. FUENTE 3 - OneWeb: validacion cruzada de arquitectura y de topologia
%  ========================================================================
fprintf('\n=======================================================================\n');
fprintf('  3. ONEWEB | hexagonal (19 haces) vs reticula 4x4 (16 haces)\n');
fprintf('=======================================================================\n');

OW = load(fullfile(ROOTDIR,'oneweb_results.mat'), 'OUT','OUTR','nm','CASES','cfg','cfgR', ...
          'RE','RD','ELIST','RLIST');

% Se emiten TODAS las entradas de crit_ratio_ffr que no son NaN, es decir los tres
% esquemas 'ffr' de CASES. En reuse1/reuseD el motor ya devuelve NaN a proposito
% (fuera del dominio de validez del criterio: no tienen sub-banda interior).
owArq = {'OneWeb-hex', 'OneWeb-rect'};
owOut = {OW.OUT, OW.OUTR};

fprintf('\n%-13s %-14s %8s %10s %12s %12s %14s %9s\n', ...
    'arquitectura','esquema','nBeams','elev med','SE_centro','SE_borde','D*SEc/SEb','renta?');
fprintf('%s\n', repmat('-',1,100));

critOW = cell(1,2);
for a = 1:2
    arq = owArq{a};  O = owOut{a};
    critOW{a} = O.diag.crit_ratio_ffr;
    for c = 1:numel(OW.nm)
        cr = O.diag.crit_ratio_ffr(c);
        if isnan(cr), continue; end               % esquema sin sub-banda interior
        SEc = O.K{c}.center.SE_mean_bpsHz;
        SEb = O.K{c}.edge.SE_mean_bpsHz;

        ROW = addrow(ROW, 'oneweb', arq, OW.nm{c}, NaN, 'nominal', O.cases{c}.Delta, O.T, ...
                     O.nBeams, O.diag.elev_mean, O.K{c}.all.n, SEc, SEb, cr);

        fprintf('%-13s %-14s %8d %10.2f %12.6f %12.6f %14.6f %9s\n', ...
            arq, OW.nm{c}, O.nBeams, O.diag.elev_mean, SEc, SEb, cr, rentaTxt(cr));
    end
end

%% --- 4b. OneWeb-hex | barrido de MASCARA DE ELEVACION (RE, 6 puntos) ----
%  RE recorre cfg.geom.minElev = 70..45 deg con el layout NOMINAL (19 haces), asi
%  que T y nBeams son los de OUT y lo unico que cambia es la geometria de servicio.
fprintf('\n--- 3b. OneWeb-hex | barrido de MASCARA DE ELEVACION (RE, %d puntos) ---\n', ...
    numel(OW.ELIST));
fprintf(['El punto de 55 deg es la mascara NOMINAL y por tanto reproduce OUT; se\n' ...
         'conserva porque el eje se lee como serie (se comprueba en la Sec. 5f).\n\n']);
fprintf('%-9s %-14s %8s %10s %12s %12s %14s %9s\n', ...
    'mascara','esquema','nBeams','elev med','SE_centro','SE_borde','D*SEc/SEb','renta?');
fprintf('%s\n', repmat('-',1,96));

critRE = nan(numel(OW.ELIST), numel(OW.nm));
for i = 1:numel(OW.ELIST)
    O = OW.RE{i};
    critRE(i,:) = O.diag.crit_ratio_ffr(:).';
    for c = 1:numel(OW.nm)
        cr = O.diag.crit_ratio_ffr(c);
        if isnan(cr), continue; end               % esquema sin sub-banda interior
        SEc = O.K{c}.center.SE_mean_bpsHz;
        SEb = O.K{c}.edge.SE_mean_bpsHz;

        ROW = addrow(ROW, 'oneweb-RE', 'OneWeb-hex', ...
                     sprintf('%s | mascara %d', OW.nm{c}, OW.ELIST(i)), ...
                     OW.ELIST(i), sprintf('mascara %d', OW.ELIST(i)), ...
                     O.cases{c}.Delta, O.T, O.nBeams, O.diag.elev_mean, ...
                     O.K{c}.all.n, SEc, SEb, cr);

        fprintf('%-9d %-14s %8d %10.2f %12.6f %12.6f %14.6f %9s\n', ...
            OW.ELIST(i), OW.nm{c}, O.nBeams, O.diag.elev_mean, SEc, SEb, cr, rentaTxt(cr));
    end
end

%% --- 4c. OneWeb-hex | barrido de ANILLOS con la HUELLA FIJA (RD) --------
%  RD densifica DENTRO de la misma huella (s = 570/R km), luego nBeams recorre
%  19/37/61/91 y el ancho de haz se estrecha con R.
%
%  RD(1) SE SALTA: nRings=2 ES la configuracion nominal y duplica exactamente la
%  fila de OUT. El assert de abajo lo verifica en vez de darlo por supuesto -- es
%  una regresion GRATIS entre dos caminos de codigo del mismo runner (OUT sale de
%  la cfg nominal; RD(1) de un cfg_over que recalcula beamwidth3dB_deg desde la
%  huella y vuelve a caer en el mismo sitio). Si algun dia dejara de cumplirse,
%  significaria que el cfg_over de RD ya no reproduce el nominal, y eso hay que
%  saberlo ANTES de emitir la fila, no despues.
crRD1 = OW.RD{1}.diag.crit_ratio_ffr(:);
crOUT = OW.OUT.diag.crit_ratio_ffr(:);
difRD1 = max(abs(crRD1(~isnan(crRD1)) - crOUT(~isnan(crOUT))));
CHK.dif_RD1_vs_OUT = difRD1;
if ~isequal(isnan(crRD1(:)), isnan(crOUT(:))) || difRD1 > 1e-12
    error('export_h2_criterion:RD1', ...
        ['RD(1) (nRings=2) deberia duplicar exactamente la fila nominal de OUT y ' ...
         'difiere en max %.3e (> 1e-12): el cfg_over del barrido de anillos ya no ' ...
         'reproduce la configuracion nominal.'], difRD1);
end
fprintf('\n--- 3c. OneWeb-hex | barrido de ANILLOS con la HUELLA FIJA (RD) ---\n');
fprintf(['[assert] RD(1) (nRings=2) == OUT (nominal): max|dif| = %.3e <= 1e-12  OK\n' ...
         '         -> RD(1) se SALTA para no duplicar la fila nominal.\n\n'], difRD1);
fprintf('%-9s %-14s %8s %10s %12s %12s %14s %9s\n', ...
    'nRings','esquema','nBeams','elev med','SE_centro','SE_borde','D*SEc/SEb','renta?');
fprintf('%s\n', repmat('-',1,96));

critRD = nan(numel(OW.RLIST), numel(OW.nm));
for i = 1:numel(OW.RLIST)
    O = OW.RD{i};
    critRD(i,:) = O.diag.crit_ratio_ffr(:).';
    if i == 1, continue; end                      % nominal: ya emitida como OUT
    for c = 1:numel(OW.nm)
        cr = O.diag.crit_ratio_ffr(c);
        if isnan(cr), continue; end
        SEc = O.K{c}.center.SE_mean_bpsHz;
        SEb = O.K{c}.edge.SE_mean_bpsHz;

        ROW = addrow(ROW, 'oneweb-RD', 'OneWeb-hex', ...
                     sprintf('%s | nRings=%d', OW.nm{c}, OW.RLIST(i)), ...
                     OW.RLIST(i), sprintf('nRings=%d', OW.RLIST(i)), ...
                     O.cases{c}.Delta, O.T, O.nBeams, O.diag.elev_mean, ...
                     O.K{c}.all.n, SEc, SEb, cr);

        fprintf('%-9d %-14s %8d %10.2f %12.6f %12.6f %14.6f %9s\n', ...
            OW.RLIST(i), OW.nm{c}, O.nBeams, O.diag.elev_mean, SEc, SEb, cr, rentaTxt(cr));
    end
end

% DIAGNOSTICO DEL EJE DE ANILLOS: por que truncar el cluster DEPRIME el cociente.
% Se lee de RD(1) y RD(end), es decir de 2 a 5 anillos con la huella FIJA.
iC3 = 4;                                          % ffr(3,a*) -> el centro es comun
if ~strcmpi(OW.CASES{iC3}.scheme,'ffr') || OW.CASES{iC3}.Delta ~= 3
    error('export_h2_criterion:iC3', ...
        'El caso %d de OneWeb es "%s" (scheme=%s, Delta=%d) y se esperaba ffr con Delta=3.', ...
        iC3, OW.nm{iC3}, OW.CASES{iC3}.scheme, OW.CASES{iC3}.Delta);
end
SEc2 = OW.RD{1}.K{iC3}.center.SE_mean_bpsHz;   SEc5 = OW.RD{end}.K{iC3}.center.SE_mean_bpsHz;
SEb2 = OW.RD{1}.K{iC3}.edge.SE_mean_bpsHz;     SEb5 = OW.RD{end}.K{iC3}.edge.SE_mean_bpsHz;
CHK.RD_SEc_2a5 = [SEc2 SEc5];
CHK.RD_SEb_2a5 = [SEb2 SEb5];
fprintf(['\n[diagnostico] De %d a %d anillos con la huella FIJA (esquema %s):\n' ...
         '   SE_centro %.5f -> %.5f  (%+.2f %%)\n' ...
         '   SE_borde  %.5f -> %.5f  (%+.2f %%)\n' ...
         '   El centro apenas se mueve y el BORDE es el que se hunde, luego anadir\n' ...
         '   anillos SUBE el cociente por el DENOMINADOR. Dicho al reves: truncar el\n' ...
         '   cluster DEPRIME crit en vez de inflarlo, que es lo contrario de lo que\n' ...
         '   se sospechaba. Importa para leer la Sec. 6: los cruces de RD NO son un\n' ...
         '   artefacto de haber recortado el cluster.\n'], ...
    OW.RLIST(1), OW.RLIST(end), OW.nm{iC3}, ...
    SEc2, SEc5, 100*(SEc5/SEc2 - 1), SEb2, SEb5, 100*(SEb5/SEb2 - 1));

%% ========================================================================
%  5. TABLA PLANA
%  ========================================================================
cruza_uno = ROW.crit > 1;                         % NaN -> false (no evaluable)
H2 = table(string(ROW.fuente(:)), string(ROW.arquitectura(:)), string(ROW.caso(:)), ...
           ROW.variante(:), ROW.Delta(:), ROW.T(:), ROW.nBeams(:), ROW.elev_mean(:), ...
           ROW.n_muestras(:), ROW.SE_c(:), ROW.SE_b(:), ROW.crit(:), cruza_uno(:), ...
    'VariableNames', {'fuente','arquitectura','caso','variante','Delta','T','nBeams', ...
                      'elev_mean','n_muestras','SE_c','SE_b','crit','cruza_uno'});

fprintf('\n=======================================================================\n');
fprintf('  4. TABLA CONSOLIDADA (%d filas)\n', height(H2));
fprintf('=======================================================================\n');
fprintf(['NOTA sobre las columnas: a las diez pedidas se anaden "caso", "variante" y\n' ...
         '"n_muestras". "caso" es IMPRESCINDIBLE (sin el, los tres tramos de E5, las tres\n' ...
         'columnas de cada T del barrido y los 6+3 puntos de RE/RD serian filas\n' ...
         'indistinguibles); "variante" lleva el valor del eje barrido -- mascara en GRADOS\n' ...
         'en las filas oneweb-RE y numero de ANILLOS en las oneweb-RD, NaN en el resto\n' ...
         '(las dos unidades comparten columna porque son ejes mutuamente excluyentes y\n' ...
         '"caso" dice siempre cual es); "n_muestras" documenta el tamano de cada\n' ...
         'submuestra, en particular el tramo VACIO de E5.\n' ...
         '"cruza_uno" es false tambien cuando crit es NaN (no evaluable, no "no cruza").\n\n']);
fprintf('%-10s %-13s %-26s %8s %6s %6s %7s %9s %10s %10s %10s %10s %5s\n', ...
    'fuente','arquitect.','caso','variante','Delta','T','nBeams','elev_med','n_muestr', ...
    'SE_c','SE_b','crit','>1?');
fprintf('%s\n', repmat('-',1,142));
for r = 1:height(H2)
    fprintf('%-10s %-13s %-26s %8s %6d %6d %7d %9s %10d %10s %10s %10s %5s\n', ...
        H2.fuente(r), H2.arquitectura(r), H2.caso(r), fmtnum(H2.variante(r),'%.0f'), ...
        H2.Delta(r), H2.T(r), H2.nBeams(r), ...
        fmtnum(H2.elev_mean(r),'%.2f'), H2.n_muestras(r), ...
        fmtnum(H2.SE_c(r),'%.5f'), fmtnum(H2.SE_b(r),'%.5f'), ...
        fmtnum(H2.crit(r),'%.6f'), ternstr(H2.cruza_uno(r),'SI',''));
end
fprintf('%s\n', repmat('-',1,142));

%% ========================================================================
%  6. VERIFICACIONES BLOQUEANTES (dos tolerancias, declaradas por bloque)
%  ========================================================================
fprintf('\n=======================================================================\n');
fprintf('  5. VERIFICACION CONTRA LAS CIFRAS PUBLICADAS\n');
fprintf('=======================================================================\n');
fprintf(['QUE TOLERANCIA USA CADA BLOQUE (y por que no es la misma):\n' ...
         '  (a) E5 ratio_bin/ratio_all ........ %.0e  referencia heredada a 6 decimales\n' ...
         '  (b) barrido de densidad, D=3 y D=4  %.0e  referencia a PRECISION COMPLETA\n' ...
         '  (c) OneWeb OUT / OUTR ............. %.0e  referencia a PRECISION COMPLETA\n' ...
         '  (f) OneWeb RE (mascara) ........... %.0e  referencia a PRECISION COMPLETA\n' ...
         '  (g) OneWeb RD (anillos) ........... %.0e  referencia a PRECISION COMPLETA\n' ...
         '  (d) y (e) controles estructurales . %.0e  el valor medido es ~1e-16\n' ...
         'Una referencia escrita a 6 decimales NO puede comprobarse a %.0e: su propio\n' ...
         'redondeo vale ya 5e-7, luego el margen de la puerta seria del orden del\n' ...
         'redondeo de la propia referencia y la comprobacion no informaria de nada.\n'], ...
         TOL, TOL9, TOL9, TOL9, TOL9, TOL, TOL9);
nChecks = 0;

% --- (a) E5 -------------------------------------------------------------
fprintf('\n(a) E5 | ratio_bin y ratio_all   [tolerancia %.0e]\n', TOL);
nChecks = nChecks + 1;
if ~isnan(E5.ratio_bin(1))
    fail{end+1} = sprintf('E5 ratio_bin(25-35 deg) = %.6f y se esperaba NaN (tramo vacio)', ...
        E5.ratio_bin(1));
    fprintf('    %-34s %12s %12s   %s\n', 'ratio_bin(25-35 deg)', 'NaN', ...
        sprintf('%.6f',E5.ratio_bin(1)), '*** FALLA ***');
else
    fprintf('    %-34s %12s %12s   %s\n', 'ratio_bin(25-35 deg)', 'NaN', 'NaN', 'OK');
end
CHK.dif_E5 = nan(1,3);
[fail, CHK.dif_E5(1)] = chk(fail, 'E5 ratio_bin(35-60 deg)', E5.ratio_bin(2), 0.777904, TOL);
[fail, CHK.dif_E5(2)] = chk(fail, 'E5 ratio_bin(60-90 deg)', E5.ratio_bin(3), 0.862729, TOL);
[fail, CHK.dif_E5(3)] = chk(fail, 'E5 ratio_all',            E5.ratio_all,    0.855612, TOL);
nChecks = nChecks + 3;

% --- (b) barrido de densidad, columnas Delta=3 y Delta=4 ----------------
%  REFERENCIAS A PRECISION COMPLETA (leidas del .mat), no a 4 decimales. La version
%  anterior comparaba [0.8461 ... 1.0202] contra 1e-4 y las desviaciones medidas
%  (1e-6 a 3.9e-5) eran del orden del propio redondeo de la referencia: la puerta
%  no distinguia una regresion real de un artefacto de transcripcion. Con la
%  referencia completa la tolerancia baja cinco ordenes, a 1e-9.
fprintf('\n(b) Barrido de densidad | columnas FFR(D=3) y FFR(D=4), %d puntos   [tolerancia %.0e]\n', ...
    nPt, TOL9);
refD4 = [0.8461390958323614
         0.8314986257192505
         0.8433668553075089
         0.8828967569397622
         0.9492638724425534
         0.9967895517833226
         1.0201783066999260];
refD3 = [0.8847155290581836
         0.8738992845516151
         0.8529940001154039
         0.8374053190915287
         0.8556122026811039
         0.8772364448650977
         0.8906985238622193];
if nPt ~= numel(refD4) || nPt ~= numel(refD3)
    error('export_h2_criterion:nPt', ...
        'El barrido tiene %d puntos y las referencias %d (D=4) / %d (D=3).', ...
        nPt, numel(refD4), numel(refD3));
end
CHK.dif_dens_D3 = nan(1,nPt);
CHK.dif_dens_D4 = nan(1,nPt);
for i = 1:nPt
    [fail, CHK.dif_dens_D3(i)] = chk(fail, sprintf('crit D=3 @ T=%d', DS.Tlist(i)), ...
        critDens(i,1), refD3(i), TOL9);
end
for i = 1:nPt
    [fail, CHK.dif_dens_D4(i)] = chk(fail, sprintf('crit D=4 @ T=%d', DS.Tlist(i)), ...
        critDens(i,2), refD4(i), TOL9);
end
nChecks = nChecks + 2*nPt;

% --- (c) OneWeb OUT / OUTR ----------------------------------------------
fprintf('\n(c) OneWeb | hexagonal (19 haces) y reticula 4x4 (16 haces)   [tolerancia %.0e]\n', TOL9);
CHK.dif_ow = nan(1,4);
[fail, CHK.dif_ow(1)] = chk(fail, 'OneWeb hex  ffr(D=3)',  critOW{1}(4), 0.8004210191013170, TOL9);
[fail, CHK.dif_ow(2)] = chk(fail, 'OneWeb hex  ffr(D=4)',  critOW{1}(5), 0.9296378604696611, TOL9);
[fail, CHK.dif_ow(3)] = chk(fail, 'OneWeb rect ffr(D=3)',  critOW{2}(4), 1.2185074307527284, TOL9);
[fail, CHK.dif_ow(4)] = chk(fail, 'OneWeb rect ffr(D=4)',  critOW{2}(5), 1.1707428939433182, TOL9);
nChecks = nChecks + 4;

% --- (d) EL CRITERIO NO DEPENDE DE LA POLITICA --------------------------
%  La columna 4 (FFR estatica Delta=3) y la 6 (FFR ADAPTATIVA Delta=3) deben dar
%  el MISMO cociente en los 7 puntos. No es una casualidad numerica ni una
%  comprobacion de rutina: es una propiedad ESTRUCTURAL que conviene poder
%  enunciar en la defensa. Tres eslabones, todos verificables en el motor:
%    (1) SE_mean_bpsHz = mean(log2(1+SINR)) (compute_kpis): eficiencia espectral
%        POR HERCIO, sin B_user dentro de la formula;
%    (2) la SINR es a su vez invariante al ancho asignado, porque la convencion es
%        PSD CONSTANTE (ffr_context) y entonces C, N e I_intra escalan TODOS con
%        B_user y el ancho se cancela en el cociente. Luego SE no depende de alpha;
%    (3) la clasificacion centro/borde la fija tau (cuantil sobre la SINR de
%        referencia en reuso-1), no alpha, asi que las dos poblaciones sobre las
%        que se promedia son las MISMAS en la estatica y en la adaptativa.
%  Consecuencia de defensa: adaptar alpha no puede cambiar el criterio que decide
%  si adaptar alpha sirve de algo. El diagnostico de H2 es PREVIO a la politica.
%  MATIZ NUMERICO: la cancelacion del eslabon (2) es EXACTA en el algebra pero
%  pasa por conversiones a dB (log10 y 10^), asi que las dos columnas coinciden a
%  PRECISION DE MAQUINA (~1e-16), no bit a bit. Por eso el criterio de aceptacion
%  es la tolerancia, y la identidad bit a bit se reporta aparte como informacion.
fprintf('\n(d) EL CRITERIO NO DEPENDE DE LA POLITICA: columna FFR(D=3) vs ADAPTATIVA(D=3)\n');
difPol = abs(critDens(:,1) - critDens(:,3));
CHK.max_dif_estatica_vs_adaptativa = max(difPol);
CHK.identicas_estatica_adaptativa  = isequaln(critDens(:,1), critDens(:,3));
fprintf('    %-6s %16s %16s %14s\n', 'T', 'FFR(D=3)', 'ADAPT(D=3)', '|dif|');
for i = 1:nPt
    fprintf('    %-6d %16.9f %16.9f %14.3e\n', ...
        DS.Tlist(i), critDens(i,1), critDens(i,3), difPol(i));
end
fprintf('    max|dif| = %.3e  (%d de %d puntos coinciden bit a bit; isequaln del vector: %s)\n', ...
    CHK.max_dif_estatica_vs_adaptativa, sum(difPol == 0), nPt, ...
    ternstr(CHK.identicas_estatica_adaptativa,'SI','NO'));
nChecks = nChecks + 1;
if CHK.max_dif_estatica_vs_adaptativa > TOL
    fail{end+1} = sprintf(['la columna FFR(D=3) y la ADAPTATIVA(D=3) del barrido difieren ' ...
        'en max %.3e (> %.0e)'], CHK.max_dif_estatica_vs_adaptativa, TOL);
else
    fprintf(['    LECTURA: las dos columnas COINCIDEN en los %d puntos (max|dif| = %.1e, es\n' ...
             '    decir precision de maquina; NO bit a bit, ver el matiz numerico del codigo).\n' ...
             '    El criterio se define sobre SE (eficiencia espectral POR HERCIO, y la SINR\n' ...
             '    es invariante al ancho por la convencion de PSD constante) y sobre la\n' ...
             '    clasificacion de REFERENCIA (fijada por tau, no por alpha), luego es\n' ...
             '    INDEPENDIENTE de la politica: la FFR adaptativa no puede mover el cociente\n' ...
             '    que decide si adaptar renta. La refutacion de H2 es PREVIA a la regla.\n'], ...
             nPt, CHK.max_dif_estatica_vs_adaptativa);
end

% --- (e) ancla cruzada E5 <-> barrido en el punto comun T=1584 -----------
%  E5 y el barrido evaluan el MISMO escenario (T=1584/72, misma rejilla y misma
%  ventana) por caminos de codigo DISTINTOS (E5 via build_ctx_blocked con
%  timeBlock=10; el barrido via run_one_density con timeBlock=4) y con alpha
%  distinto (alpha* de E5 frente a 0.4). Que los dos SE coincidan es a la vez una
%  regresion entre runners y una confirmacion empirica del punto (d).
iT1584 = find(DS.Tlist == T_E5, 1);
CHK.ancla_E5_vs_densidad = NaN;
if ~isempty(iT1584)
    SEc_ds = DS.results{iT1584}.K{4}.center.SE_mean_bpsHz;
    SEb_ds = DS.results{iT1584}.K{4}.edge.SE_mean_bpsHz;
    CHK.ancla_E5_vs_densidad = max(abs([SEc_ds - SEc_E5, SEb_ds - SEb_E5, ...
                                        critDens(iT1584,1) - E5.ratio_all]));
    fprintf('\n(e) ANCLA CRUZADA E5 <-> barrido en el punto comun T=%d\n', T_E5);
    fprintf('    %-12s %18s %18s %14s\n', 'magnitud', 'E5', 'barrido', '|dif|');
    fprintf('    %-12s %18.9f %18.9f %14.3e\n', 'SE_centro', SEc_E5, SEc_ds, abs(SEc_ds-SEc_E5));
    fprintf('    %-12s %18.9f %18.9f %14.3e\n', 'SE_borde',  SEb_E5, SEb_ds, abs(SEb_ds-SEb_E5));
    fprintf('    %-12s %18.9f %18.9f %14.3e\n', 'crit', E5.ratio_all, critDens(iT1584,1), ...
        abs(critDens(iT1584,1)-E5.ratio_all));
    fprintf('    max|dif| = %.3e  (con alpha DISTINTO: %.2f en E5 frente a %.2f en el barrido)\n', ...
        CHK.ancla_E5_vs_densidad, E5.cases{iST}.alpha, DS.schemeCases{4}.alpha);
    nChecks = nChecks + 1;
    if CHK.ancla_E5_vs_densidad > TOL
        fail{end+1} = sprintf(['E5 y el barrido discrepan en el punto comun T=%d ' ...
            '(max %.3e)'], T_E5, CHK.ancla_E5_vs_densidad);
    end
end

% --- (f) OneWeb RE | barrido de mascara de elevacion --------------------
fprintf('\n(f) OneWeb-hex | barrido de MASCARA (RE), ffr(D=3) y ffr(D=4)   [tolerancia %.0e]\n', TOL9);
refRE = [70  0.8292827151560956  1.0179809043208412
         65  0.8096031022645123  0.9678257063313352
         60  0.8008990954733524  0.9372055075591955
         55  0.8004210191013170  0.9296378604696611
         50  0.8006270672125827  0.9295587535260347
         45  0.8006270840624649  0.9295587718552042];
if ~isequal(OW.ELIST(:), refRE(:,1))
    error('export_h2_criterion:ELIST', ...
        'ELIST del .mat es [%s] y la referencia [%s].', ...
        num2str(OW.ELIST(:).'), num2str(refRE(:,1).'));
end
CHK.dif_RE = nan(size(refRE,1), 2);
for i = 1:size(refRE,1)
    [fail, CHK.dif_RE(i,1)] = chk(fail, sprintf('RE ffr(D=3) | mascara %d', refRE(i,1)), ...
        critRE(i,4), refRE(i,2), TOL9);
    [fail, CHK.dif_RE(i,2)] = chk(fail, sprintf('RE ffr(D=4) | mascara %d', refRE(i,1)), ...
        critRE(i,5), refRE(i,3), TOL9);
end
nChecks = nChecks + 2*size(refRE,1);

% Cada punto lleva su PROPIA referencia y por eso el aviso siguiente NO se ha
% convertido en una comprobacion de igualdad entre puntos: las mascaras de 50 y 45
% deg dan KPIs identicos (por debajo de 55 la mascara deja de morder) pero su crit
% difiere en ~1.7e-8. Es ruido numerico del pipeline, no un fallo, y exigir
% igualdad EXACTA entre esas dos filas convertiria ese ruido en una falsa alarma.
CHK.dif_RE_50_vs_45 = max(abs(critRE(5,[4 5]) - critRE(6,[4 5])));
fprintf(['    [aviso] mascaras 50 y 45 deg: KPIs identicos pero crit difiere en\n' ...
         '            %.2e (ruido numerico, NO un fallo). No se comprueba como\n' ...
         '            igualdad entre puntos; cada uno lleva su propia referencia.\n'], ...
    CHK.dif_RE_50_vs_45);

% CONTROL EXTRA (regresion gratis): el punto de 55 deg de RE es la mascara NOMINAL
% y debe reproducir OUT. Se conserva en la tabla, a diferencia de RD(1), porque el
% eje de mascara se lee como serie y quitarle un punto interior lo mutilaria.
iN55 = find(OW.ELIST == OW.cfg.geom.minElev, 1);
CHK.dif_RE_nominal_vs_OUT = NaN;
if ~isempty(iN55)
    cRE = critRE(iN55,:).';
    CHK.dif_RE_nominal_vs_OUT = max(abs(cRE(~isnan(cRE)) - crOUT(~isnan(crOUT))));
    [fail, ~] = chkval(fail, sprintf('RE @ mascara nominal %d == OUT', OW.cfg.geom.minElev), ...
        CHK.dif_RE_nominal_vs_OUT, TOL9);
    nChecks = nChecks + 1;
end

% --- (g) OneWeb RD | barrido de anillos con la huella fija --------------
fprintf('\n(g) OneWeb-hex | barrido de ANILLOS (RD), ffr(D=3) y ffr(D=4)   [tolerancia %.0e]\n', TOL9);
refRD = [2  19  0.8004210191013170  0.9296378604696611
         3  37  0.9141187677667197  1.0917039841607319
         4  61  0.9212508993937801  1.0885155297881888
         5  91  0.9417007306872324  1.1052016728220038];
if ~isequal(OW.RLIST(:), refRD(:,1))
    error('export_h2_criterion:RLIST', ...
        'RLIST del .mat es [%s] y la referencia [%s].', ...
        num2str(OW.RLIST(:).'), num2str(refRD(:,1).'));
end
CHK.dif_RD = nan(size(refRD,1), 2);
for i = 1:size(refRD,1)
    % nBeams = 3R(R+1)+1: se comprueba de paso que el cfg_over de anillos hizo lo
    % que dice (si nRings no se hubiera aplicado, los 4 puntos tendrian 19 haces).
    if OW.RD{i}.nBeams ~= refRD(i,2)
        fail{end+1} = sprintf('RD nRings=%d tiene %d haces y se esperaban %d', ...
            refRD(i,1), OW.RD{i}.nBeams, refRD(i,2));
    end
    sfx = ternstr(i == 1, ' [nominal]', '');
    [fail, CHK.dif_RD(i,1)] = chk(fail, sprintf('RD ffr(D=3) | nRings=%d%s', refRD(i,1), sfx), ...
        critRD(i,4), refRD(i,3), TOL9);
    [fail, CHK.dif_RD(i,2)] = chk(fail, sprintf('RD ffr(D=4) | nRings=%d%s', refRD(i,1), sfx), ...
        critRD(i,5), refRD(i,4), TOL9);
end
nChecks = nChecks + 2*size(refRD,1) + size(refRD,1);

% --- Veredicto de la verificacion ---------------------------------------
CHK.tol      = TOL;
CHK.tol_crit = TOL9;
CHK.nChecks  = nChecks;
CHK.nFail    = numel(fail);
fprintf('\n-----------------------------------------------------------------------\n');
if ~isempty(fail)
    fprintf('VERIFICACION FALLIDA: %d de %d comprobaciones fuera de tolerancia.\n', ...
        numel(fail), nChecks);
    fprintf('  - %s\n', fail{:});
    fprintf('NO se guarda h2_criterion.mat.\n');
    error('export_h2_criterion:verificacion', ...
        ['%d comprobacion(es) no coinciden con las cifras publicadas. Ver el log %s.'], ...
        numel(fail), LOGF);
end
fprintf(['VERIFICACION SUPERADA: las %d cifras de control coinciden (crit a %.0e con\n' ...
         'referencia de precision completa; referencias heredadas y controles\n' ...
         'estructurales a %.0e).\n'], nChecks, TOL9, TOL);
fprintf('-----------------------------------------------------------------------\n');

%% ========================================================================
%  7. LO QUE CRUZA 1 (el hallazgo)
%  ========================================================================
fprintf('\n=======================================================================\n');
fprintf('  6. CONFIGURACIONES CON  Delta*SE_centro/SE_borde > 1  (el interior RENTA)\n');
fprintf('=======================================================================\n');

iX = find(H2.cruza_uno);
fprintf('Filas de la tabla con crit > 1: %d de %d.\n\n', numel(iX), height(H2));
for k = 1:numel(iX)
    r = iX(k);
    fprintf('  [%d] %-13s %-26s Delta=%d  T=%-5d nBeams=%2d  ->  crit = %.6f\n', ...
        k, H2.arquitectura(r), H2.caso(r), H2.Delta(r), H2.T(r), H2.nBeams(r), H2.crit(r));
end

% CONFIGURACIONES DISTINTAS. La deduplicacion existe porque la FFR adaptativa
% duplica el numero de su estatica del mismo Delta (punto (d)) y no es un
% contraejemplo independiente: esas dos SI deben colapsar.
%
% LA CLAVE ANTIGUA ERA (arquitectura, T, Delta) Y SE QUEDABA CORTA al entrar RE y
% RD: todas las filas de OneWeb-hex comparten arquitectura y T=720, luego las de
% Delta=4 que cruzan (mascara 70, nRings 3, 4 y 5) colapsaban en UNA sola y el
% resumen informaba de un cruce donde hay cuatro configuraciones fisicas distintas.
% La clave lleva ahora tambien nBeams y el PUNTO del eje ('eje': la parte de 'caso'
% sin la politica), que es lo unico que separa la mascara de 70 deg de la nominal
% -- las dos tienen 19 haces y el mismo T. Al excluir la politica del identificador,
% estatica y adaptativa del mismo Delta siguen colapsando, que es lo que se queria.
ejeAll = string(ROW.eje(:));
keyAll = strcat(H2.arquitectura, "|T=", string(H2.T), "|D=", string(H2.Delta), ...
                "|nB=", string(H2.nBeams), "|", ejeAll);
CHK.nFilasCruzan  = numel(iX);
CHK.nConfigCruzan = 0;
if ~isempty(iX)
    [uk, ia] = unique(keyAll(iX), 'stable');
    CHK.nConfigCruzan = numel(uk);
    fprintf(['\nCONFIGURACIONES DISTINTAS que cruzan 1: %d\n' ...
             '(clave = arquitectura | T | Delta | nBeams | punto del eje, SIN la politica)\n'], ...
        numel(uk));
    for k = 1:numel(uk)
        r = iX(ia(k));
        fprintf('  * %-13s T=%-5d Delta=%d  nBeams=%2d  %-14s  SE_c=%.5f  SE_b=%.5f  crit = %.6f\n', ...
            H2.arquitectura(r), H2.T(r), H2.Delta(r), H2.nBeams(r), ejeAll(r), ...
            H2.SE_c(r), H2.SE_b(r), H2.crit(r));
    end
    nDup = numel(iX) - numel(uk);
    if nDup > 0
        fprintf(['  (La%s %d fila%s restante%s %s la FFR ADAPTATIVA del mismo Delta, que da el\n' ...
                 '   MISMO cociente por el punto (d): no %s contraejemplo%s independiente%s.)\n'], ...
            ternstr(nDup==1,'',  's'), nDup, ternstr(nDup==1,'','s'), ternstr(nDup==1,'','s'), ...
            ternstr(nDup==1,'es','son'), ternstr(nDup==1,'es','son'), ...
            ternstr(nDup==1,'','s'), ternstr(nDup==1,'','s'));
    end
end

%% --- RECUENTO DESAGREGADO POR (Delta, topologia) ------------------------
%  El enunciado que va a la memoria es un patron, no una lista de filas, asi que el
%  log tiene que sostenerlo sin que nadie cuente a mano. Se cuenta sobre
%  CONFIGURACIONES DEDUPLICADAS y EVALUABLES (crit no NaN): una fila con crit = NaN
%  no es "no cruza", es "no se puede decidir", y mezclarlas falsearia el
%  denominador. La topologia se deriva de la arquitectura: hex / rect.
topo = repmat("hex", height(H2), 1);
topo(contains(H2.arquitectura, "rect")) = "rect";
[~, iU]  = unique(keyAll, 'stable');              % una fila representante por config
okU      = iU(~isnan(H2.crit(iU)));               % configuraciones EVALUABLES
fprintf('\n-----------------------------------------------------------------------\n');
fprintf('RECUENTO DE CRUCES POR (Delta, topologia)\n');
fprintf('[configuraciones deduplicadas y EVALUABLES; una fila con crit NaN no cuenta]\n\n');
fprintf('%-10s %7s %13s %9s   %s\n', 'topologia','Delta','evaluables','cruzan','arquitecturas');
fprintf('%s\n', repmat('-',1,78));
tl = ["hex" "rect"];  dl = unique(H2.Delta(okU)).';
cT = strings(0,1);  cD = [];  cN = [];  cX = [];  cA = strings(0,1);
for tt = tl
    for dd = dl
        sel = okU(topo(okU) == tt & H2.Delta(okU) == dd);
        if isempty(sel), continue; end
        nCr  = sum(H2.cruza_uno(sel));
        arqs = strjoin(unique(H2.arquitectura(sel), 'stable').', ', ');
        fprintf('%-10s %7d %13d %9d   %s\n', tt, dd, numel(sel), nCr, arqs);
        cT(end+1,1) = tt;   cD(end+1,1) = dd;   %#ok<AGROW>
        cN(end+1,1) = numel(sel);  cX(end+1,1) = nCr;  cA(end+1,1) = string(arqs); %#ok<AGROW>
    end
end
fprintf('%s\n', repmat('-',1,78));
CNT = table(cT, cD, cN, cX, cA, 'VariableNames', ...
    {'topologia','Delta','evaluables','cruzan','arquitecturas'});

nH3  = sum(CNT.evaluables(CNT.topologia=="hex"  & CNT.Delta==3));
xH3  = sum(CNT.cruzan(    CNT.topologia=="hex"  & CNT.Delta==3));
nH4  = sum(CNT.evaluables(CNT.topologia=="hex"  & CNT.Delta==4));
xH4  = sum(CNT.cruzan(    CNT.topologia=="hex"  & CNT.Delta==4));
nR   = sum(CNT.evaluables(CNT.topologia=="rect"));
xR   = sum(CNT.cruzan(    CNT.topologia=="rect"));
CHK.cruces_hex_D3  = [xH3 nH3];
CHK.cruces_hex_D4  = [xH4 nH4];
CHK.cruces_rect    = [xR  nR];
fprintf(['\nEL ENUNCIADO QUE SOSTIENE ESTE RECUENTO:\n' ...
         '  * Delta=3 sobre reticula HEXAGONAL: el cociente NO cruza en NINGUNA de las\n' ...
         '    %d configuraciones evaluables, y eso incluye las DOS arquitecturas\n' ...
         '    (Starlink y OneWeb) y los cuatro ejes barridos (elevacion, densidad\n' ...
         '    orbital, mascara y anillos). Cruces: %d de %d.\n' ...
         '  * Delta=4 sobre reticula HEXAGONAL: cruza en %d de %d configuraciones.\n' ...
         '  * Reticula CUADRADA (4x4): cruza con LOS DOS Delta (%d de %d).\n' ...
         '  Es decir el contraejemplo fuerte no es el Delta ni la densidad: es la\n' ...
         '  TOPOLOGIA. Encaja con el hallazgo ya publicado de que la ventaja de Delta=3\n' ...
         '  es un constructo HEXAGONAL (sus co-canal caen a raiz(3)*s, y sobre celdas\n' ...
         '  cuadradas a solo raiz(2)*s).\n'], nH3, xH3, nH3, xH4, nH4, xR, nR);

fprintf(['\nLECTURA (material para la memoria):\n' ...
         '  La afirmacion "Todos < 1" de INFORME_TECNICO Sec. 9.5 NO se sostiene sobre el\n' ...
         '  corpus completo. Su tabla de evidencia cita del barrido de densidad SOLO la\n' ...
         '  columna Delta=3 y de OneWeb SOLO la topologia hexagonal NOMINAL, que son justo\n' ...
         '  los cortes en los que el cociente no llega a 1. Al mirar tambien Delta=4, la\n' ...
         '  reticula 4x4 y los dos barridos de OneWeb (mascara y anillos) aparecen los\n' ...
         '  contraejemplos de arriba.\n' ...
         '  MATIZ, para no sobreinterpretarlo: el criterio es NECESARIO, no suficiente.\n' ...
         '  Que crit > 1 dice que subir alpha aumentaria el agregado por celda, NO que la\n' ...
         '  FFR adaptativa gane; medirlo exige reejecutar E5 en esas configuraciones.\n' ...
         '  Ademas T=4000 esta muy por encima de cualquier constelacion desplegada y alli\n' ...
         '  el cruce es del 2%%, mientras que en la reticula 4x4 es del 17-22%%.\n']);

fprintf(['\n  LO QUE ANADEN LAS FILAS NUEVAS (RE y RD), y va en contra de la sospecha\n' ...
         '  inicial: al anadir anillos con la huella FIJA (de %d a %d), SE_centro se mueve\n' ...
         '  un %+.2f%% (%.5f -> %.5f) mientras SE_borde cae un %+.2f%% (%.5f -> %.5f).\n' ...
         '  El cociente sube, por tanto, por el DENOMINADOR: densificar castiga al borde\n' ...
         '  (mas haces co-canal a menor separacion fisica) y deja el interior casi igual,\n' ...
         '  porque el interior ya iba en reuso-1 y no tenia proteccion que perder.\n' ...
         '  Leido al reves, que es lo que importa: TRUNCAR el cluster DEPRIME crit en vez\n' ...
         '  de inflarlo. Los cruces de RD no son, por tanto, un artefacto de haber\n' ...
         '  recortado el layout -- que era la explicacion facil y resulta ser la contraria.\n'], ...
    OW.RLIST(1), OW.RLIST(end), 100*(SEc5/SEc2-1), SEc2, SEc5, ...
    100*(SEb5/SEb2-1), SEb2, SEb5);

%% ========================================================================
%  8. GUARDADO
%  ========================================================================
META = struct();
META.generado       = stamp;
META.script         = [mfilename '.m'];
META.fuentes        = SRC;
META.tolerancia     = TOL;
META.tolerancia_crit = TOL9;
META.criterio       = 'crit = Delta * SE_centro / SE_borde; el interior renta si crit > 1';
META.SE_definicion  = 'SE_mean_bpsHz = mean(log2(1+SINR_lin)) sobre la submuestra (compute_kpis)';
META.columnas_extra = {'caso: identificador de la fila (sin el, filas indistinguibles)', ...
                       'variante: valor del eje barrido (grados de mascara en oneweb-RE, numero de anillos en oneweb-RD, NaN en el resto)', ...
                       'n_muestras: tamano de la submuestra (usuario,instante)'};
META.dedup_clave    = 'arquitectura | T | Delta | nBeams | punto del eje (SIN la politica): la estatica y la adaptativa del mismo Delta colapsan; dos puntos distintos del mismo eje no.';
META.limitaciones   = { ...
    'SE_c y SE_b por tramo de elevacion (E5) van NaN: el .mat guarda solo ratio_bin.', ...
    'cruza_uno = false tambien cuando crit es NaN (no evaluable).', ...
    'El criterio solo aplica a esquemas ffr; el motor devuelve NaN en reuse1/reuseD.', ...
    'RD(1) (nRings=2) NO se emite: duplica la fila nominal de OUT (verificado a 1e-12).', ...
    'RE @ mascara nominal SI se emite aunque duplique OUT: el eje se lee como serie.', ...
    'variante mezcla dos unidades (grados y anillos) en una columna; caso las separa.'};

outf = fullfile(RESDIR, 'h2_criterion.mat');
save(outf, 'H2', 'CNT', 'CHK', 'META');
fprintf('\nTabla guardada en %s (%d filas x %d columnas).\n', outf, height(H2), width(H2));
fprintf('Log en %s\n', LOGF);
fprintf('\nFIN.\n');
diary off;
end

%% ========================= funciones locales =========================
function R = addrow(R, fu, ar, ca, va, ej, D, T, nB, el, n, sc, sb, cr)
%ADDROW  Anade una fila a los acumuladores de la tabla plana.
%   va = variante (valor del eje barrido, NaN si la fila no pertenece a un barrido)
%   ej = etiqueta del punto del eje SIN la politica (clave de deduplicacion)
R.fuente{end+1}       = fu;
R.arquitectura{end+1} = ar;
R.caso{end+1}         = ca;
R.variante(end+1)     = va;
R.eje{end+1}          = ej;
R.Delta(end+1)        = D;
R.T(end+1)            = T;
R.nBeams(end+1)       = nB;
R.elev_mean(end+1)    = el;
R.n_muestras(end+1)   = n;
R.SE_c(end+1)         = sc;
R.SE_b(end+1)         = sb;
R.crit(end+1)         = cr;
end

function [fail, dif] = chk(fail, nombre, valor, ref, tol)
%CHK  Compara un valor con su referencia publicada y acumula el fallo si procede.
%   Los valores se imprimen a 10 decimales (no a 6) porque parte de las referencias
%   estan transcritas a precision completa y se comprueban a 1e-9: con 6 decimales
%   el log no permitiria distinguir un fallo de un redondeo de impresion.
dif = abs(valor - ref);
ok  = dif <= tol;
fprintf('    %-38s %16.10f %16.10f  |dif| = %.2e (tol %.0e)  %s\n', ...
    nombre, ref, valor, dif, tol, ternstr(ok,'OK','*** FALLA ***'));
if ~ok
    fail{end+1} = sprintf('%s: %.12f frente a %.12f publicado (|dif| = %.2e > %.0e)', ...
        nombre, valor, ref, dif, tol);
end
end

function [fail, ok] = chkval(fail, nombre, dif, tol)
%CHKVAL  Comprueba una desviacion YA calculada (identidades entre dos lecturas del
%   mismo .mat, donde no hay "referencia publicada" sino una diferencia que debe
%   ser nula). Misma politica de acumulacion que chk.
ok = dif <= tol;
fprintf('    %-38s %16s %16.3e  %25s  %s\n', nombre, '(max|dif|)', dif, ...
    sprintf('(tol %.0e)', tol), ternstr(ok,'OK','*** FALLA ***'));
if ~ok
    fail{end+1} = sprintf('%s: max|dif| = %.3e > %.0e', nombre, dif, tol);
end
end

function s = ternstr(c, a, b)
if c, s = a; else, s = b; end
end

function s = rentaTxt(cr)
%RENTATXT  Veredicto del criterio, distinguiendo NaN (no evaluable) de "no renta".
if isnan(cr), s = 'n/d'; elseif cr > 1, s = 'SI (>1)'; else, s = 'NO'; end
end

function s = fmtnum(v, f)
%FMTNUM  Formatea un numero dejando NaN legible en las tablas de texto.
if isnan(v), s = 'NaN'; else, s = sprintf(f, v); end
end
