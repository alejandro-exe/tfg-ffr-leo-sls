%% RUN_ONEWEB_DEMO  OneWeb como SISTEMA PROPIO: validacion cruzada de arquitecturas.
%  Runner DELGADO: no contiene fisica. Todo sale del motor compartido a traves de
%  run_one_density, variando SOLO cfg (perfil config_oneweb). Es el mismo principio
%  con el que E0, run_ffr_demo, E3, E3b y E4 reutilizan la misma cadena.
%
%  PREGUNTA: las conclusiones del TFG (H1, criterio de viabilidad, compresion
%  angular, auto-similitud de la reticula) se obtuvieron todas con UN sistema
%  (Starlink-like, 550 km, 53 deg, mascara 25 deg, celda 17 km). Si se reproducen en
%  una arquitectura RADICALMENTE distinta -- OneWeb: 1200 km, 87.9 deg polar,
%  mascara 55 deg, celda 285 km -- dejan de ser propiedades del sistema simulado y
%  pasan a ser propiedades del REUSO MULTIHAZ EN LEO. Si alguna NO se reproduce, el
%  resultado acota su dominio de validez, que vale igual.
%
%  SECCIONES
%    1. Perfil, layout y guarda (dos criterios, como siempre)
%    2. AUTO-SIMILITUD: dominio de validez de la identidad u = 1.61634*sqrt(3)
%    3. E2 equivalente: 6 esquemas, KPIs y veredictos de viabilidad
%    4. E3b equivalente: barrido de mascara de elevacion
%    5. Contraste de TOPOLOGIA: reticula 4x4 (16 haces) vs hexagonal (19)
%    6. E3 equivalente: densificacion de haces con la HUELLA FIJA
%    7. Resumen comparativo Starlink vs OneWeb
%
%  Resultados en oneweb_results.mat, figuras en figs_oneweb/.
%
%  Aportacion propia (diseno del experimento de validacion cruzada).

clear; clc; close all;
if ~exist('figs_oneweb','dir'), mkdir('figs_oneweb'); end

ALPHA_STAR = 0.5;      % mismo alpha* que E3/E3b, para comparar en igualdad
TAU_Q      = 50;       % cuantil de clasificacion centro/borde (mismo que Starlink)

%% ======================= 1. PERFIL Y LAYOUT =============================
fprintf('=======================================================================\n');
fprintf('  ONEWEB COMO SISTEMA PROPIO -- validacion cruzada de arquitecturas\n');
fprintf('=======================================================================\n\n');

cfg = config_oneweb('hex');
cfg.time.dt       = 120;               % Nt = 61 (el convergido de la FASE A)
cfg.time.duration = 7200;
cfg.ffr.tau_mode  = 'quantile';
cfg.ffr.tau_q     = TAU_Q;
cfg.ffr.sched     = 'share';
cfg.compute.timeBlock = 10;            % EXACTO; acota el pico de RAM

BL    = build_beam_layout(cfg);
users = build_user_grid(cfg);
Nt    = numel(cfg.time.t0:cfg.time.dt:cfg.time.t0+cfg.time.duration);

fprintf('OPERADOR: %s | h=%d km | inc=%.1f deg | %s | T=%d P=%d | mascara %d deg\n', ...
    cfg.constellations(1).name, cfg.constellations(1).h, cfg.constellations(1).inc, ...
    cfg.constellations(1).pattern, cfg.constellations(1).T, cfg.constellations(1).P, ...
    cfg.geom.minElev);
fprintf('RADIO   : f=%.1f GHz | B=%d MHz | EIRPdens=%.2f dBW/MHz | G/T=%.1f dB/K | Grx=%.1f dBi\n', ...
    cfg.radio.freq_GHz, cfg.radio.B_MHz, cfg.radio.EIRPdensity_dBWMHz, ...
    cfg.radio.GT_dBK, cfg.radio.Grx_max_dBi);
fprintf('LAYOUT  : HPBW=%.4f deg | s=%.4f deg = %.2f km | %d haces (%d anillos)\n', ...
    cfg.radio.beamwidth3dB_deg, BL.spacing_deg, BL.spacing_km, BL.nBeams, cfg.beams.nRings);
fprintf('REJILLA : radio %d km | step %d km | M=%d usuarios | Nt=%d (dt=%d s)\n', ...
    cfg.ground.radius_km, cfg.ground.step_km, users.M, Nt, cfg.time.dt);

% --- SUPUESTO DECLARADO: anclaje en la huella, no en Gmax ---
Gimp = 10*log10(41253/cfg.radio.beamwidth3dB_deg^2);
fprintf('\n[SUPUESTO] La HPBW se ancla en la HUELLA declarada (1140 km / 4 = %.0f km de celda),\n', BL.spacing_km);
fprintf('           no en Gmax. Directividad implicita %.1f dBi vs %.1f dBi declarados (%.1f dB).\n', ...
    Gimp, cfg.radio.Gmax_dBi, cfg.radio.Gmax_dBi - Gimp);
fprintf('           Es lo correcto: la huella fija la GEOMETRIA DE LA INTERFERENCIA (lo que se\n');
fprintf('           estudia) y Gmax_dBi NO entra en la fisica (auditoria, hallazgo G2).\n');

% --- GUARDA: los dos criterios de siempre ---
dcell = sqrt(sum(BL.cell_km.^2,2));
dext  = min(dcell(BL.ring == max(BL.ring)));
guard = (dext - cfg.ground.radius_km) / BL.spacing_km;
fprintf('\n[GUARDA] anillo externo a %.0f km | rejilla %d km | guarda ANALITICA = %.2f anillos -> %s\n', ...
    dext, cfg.ground.radius_km, guard, tern(guard >= 1, 'SUFICIENTE', 'ESCASA'));
fprintf('         (la comprobacion OPERATIVA, max theta vs circunradio, la emite ffr_context)\n');

% --- POR QUE nRings = 2 ES EL MODELO CORRECTO ---
fprintf('\n[MODELO] nRings=2 (19 haces ~ los 16 declarados) NO es una aproximacion por\n');
fprintf('         conveniencia: las celdas mas alla de la huella las sirven OTROS satelites\n');
fprintf('         OneWeb, y por coordinacion INTRA-operador NO son co-canal con nosotros\n');
fprintf('         (mismo argumento que justifica satPointing=''nadir''). Subir nRings para\n');
fprintf('         "ganar guarda" INVENTARIA interferencia co-canal inexistente.\n');

% --- Coloreado ---
fprintf('\n');
for D = [3 4]
    [~, COL] = ffr_coloring(BL, D, false);
    fprintf('[coloreado] Delta=%d -> reparto %s | dmin = %.3f*s = %.0f km\n', ...
        D, mat2str(COL.count.'), COL.dmin_s, COL.dmin_km);
end

%% =============== 2. AUTO-SIMILITUD: DOMINIO DE VALIDEZ ==================
fprintf('\n=======================================================================\n');
fprintf('  2. AUTO-SIMILITUD ANGULAR: dominio de validez de la identidad\n');
fprintf('=======================================================================\n');
% MATIZ AL LEER LA SALIDA: el texto que sigue llama "vecino co-canal del 1er anillo"
% al vecino ADYACENTE (theta_off = s), que es el primer co-canal SOLO en reuso pleno.
% Con Delta=3 el primer co-canal esta a sqrt(3)*s (1.500 anchos, -17.954 dB) y con
% Delta=4 a 2*s (1.732 anchos, -18.456 dB). La invariancia se cumple igual en las tres
% distancias (recorridos 0.010 / 0.014 / 0.028 dB), luego lo que se acota es el NUMERO,
% no la propiedad. Ojo: el sqrt(3) de 1.61634*sqrt(3) es el de s = sqrt(3)*sin(HPBW/2),
% NO la distancia de reuso de Delta=3. El texto impreso no se reescribe a proposito,
% para no alterar la salida ya publicada del experimento.
fprintf(['Identidad (E3/E3b): con s = sqrt(3)*sin(HPBW/2) el vecino co-canal del 1er\n' ...
         'anillo cae en u = 1.61634*sqrt(3) = 2.79958 -> G = -10.669 dB, INDEPENDIENTE\n' ...
         'del ancho de haz. Pero esa formula de 38.821 es una aprox. de ANGULO PEQUENO,\n' ...
         'asi que la invariancia tiene un DOMINIO. OneWeb (HPBW %.1f deg) lo explora.\n\n'], ...
         cfg.radio.beamwidth3dB_deg);

HPlist = [0.5 1.0 2.08 5 10 cfg.radio.beamwidth3dB_deg 20];
Ginf   = beam_gain_dB(rad2deg(sqrt(3)*sin(deg2rad(1e-4)/2)), 1e-4);
AS = struct('HPBW',HPlist,'s_deg',nan(size(HPlist)),'u',nan(size(HPlist)), ...
            'G',nan(size(HPlist)),'dev',nan(size(HPlist)));
fprintf('%11s %11s %12s %13s %14s   %s\n','HPBW[deg]','s_deg','u(vecino)','G_vec[dB]','desv.[dB]','sistema');
for k = 1:numel(HPlist)
    sd = rad2deg(sqrt(3)*sin(deg2rad(HPlist(k))/2));
    AS.s_deg(k) = sd;
    AS.u(k)     = 1.61634*sind(sd)/sind(HPlist(k)/2);
    AS.G(k)     = beam_gain_dB(sd, HPlist(k));
    AS.dev(k)   = AS.G(k) - Ginf;
    lab = '';
    if abs(HPlist(k)-2.08) < 1e-9, lab = '<- Starlink'; end
    if abs(HPlist(k)-cfg.radio.beamwidth3dB_deg) < 1e-9, lab = '<- OneWeb'; end
    fprintf('%11.3f %11.4f %12.5f %13.3f %14.3f   %s\n', ...
        HPlist(k), sd, AS.u(k), AS.G(k), AS.dev(k), lab);
end
AS.Ginf = Ginf;  AS.uinf = 1.61634*sqrt(3);
fprintf('  asintotico (HPBW->0): u = %.5f -> G = %.3f dB\n', AS.uinf, Ginf);
fprintf(['\n[HALLAZGO] La identidad no es universal: es EXACTA en regimen de angulo\n' ...
         '  pequeno (donde s = sqrt(3)*sin(HPBW/2) es lineal) y aproximada para haz\n' ...
         '  ancho. FORMULACION RECOMENDADA: "invariante para HPBW <~ 10 deg (0.005 dB\n' ...
         '  a 2 deg), aproximada para haz ancho (0.25 dB a 15 deg)". Un resultado\n' ...
         '  ACOTADO vale mas que uno universal sin acotar.\n']);

%% ============ 3. E2 EQUIVALENTE: 6 ESQUEMAS EN EL PUNTO NOMINAL =========
fprintf('\n=======================================================================\n');
fprintf('  3. COMPARATIVA DE ESQUEMAS (mascara nominal %d deg)\n', cfg.geom.minElev);
fprintf('=======================================================================\n');

CASES = { ...
  struct('scheme','reuse1','Delta',1,'alpha',NaN,      'adaptive',false,'name','reuse1'), ...
  struct('scheme','reuseD','Delta',3,'alpha',NaN,      'adaptive',false,'name','reuseD(3)'), ...
  struct('scheme','reuseD','Delta',4,'alpha',NaN,      'adaptive',false,'name','reuseD(4)'), ...
  struct('scheme','ffr',   'Delta',3,'alpha',ALPHA_STAR,'adaptive',false,'name','ffr(3,a*)'), ...
  struct('scheme','ffr',   'Delta',4,'alpha',ALPHA_STAR,'adaptive',false,'name','ffr(4,a*)'), ...
  struct('scheme','ffr',   'Delta',3,'alpha',ALPHA_STAR,'adaptive',true, 'name','ffr-adapt(3)')};
nm = cellfun(@(c) c.name, CASES, 'UniformOutput', false);

pt = struct('cases',{CASES}, 'elMin_valid',70, 'label','OneWeb nominal', ...
            'verbose',true, 'idx',1);
tic; OUT = run_one_density(cfg, pt); tNom = toc;
fprintf('[punto nominal] %.1f s | cobertura %.3f | elev media %.1f deg | N=%d\n\n', ...
    tNom, OUT.diag.covFrac, OUT.diag.elev_mean, OUT.N);

print_scheme_table(OUT.K, nm, 'VENTANA COMPLETA (elev >= 55 deg)');
if ~isempty(OUT.KV)
    print_scheme_table(OUT.KV, nm, 'GEOMETRIA FAVORABLE (elev >= 70 deg)');
end

% H1 en el borde
r1 = OUT.K{1}; f3 = OUT.K{4}; f4 = OUT.K{5};
fprintf('\n[H1] R_p5 de BORDE: reuse1 %.2f -> ffr(D=3) %.2f Mbps (x%.2f) | ffr(D=4) %.2f (x%.2f)\n', ...
    r1.edge.R_p5_Mbps, f3.edge.R_p5_Mbps, f3.edge.R_p5_Mbps/r1.edge.R_p5_Mbps, ...
    f4.edge.R_p5_Mbps, f4.edge.R_p5_Mbps/r1.edge.R_p5_Mbps);
fprintf('[H1] SINR_edge_p5: reuse1 %+.2f -> ffr(D=3) %+.2f dB (%+.2f dB) | ffr(D=4) %+.2f (%+.2f dB)\n', ...
    r1.viab.SINR_edge_p5, f3.viab.SINR_edge_p5, f3.viab.SINR_edge_p5-r1.viab.SINR_edge_p5, ...
    f4.viab.SINR_edge_p5, f4.viab.SINR_edge_p5-r1.viab.SINR_edge_p5);
% diag.crit_ratio_ffr viene POR ESQUEMA y ya sale NaN fuera de los esquemas FFR
% (run_one_density lo restringe a su dominio de validez): el criterio
% Delta*SE_centro/SE_borde > 1 presupone una sub-banda INTERIOR, y en reuse1 y
% reuseD centro y borde ven el MISMO conjunto co-canal. El filtro isFFR se mantiene
% aqui para IMPRIMIR solo lo que aplica; la garantia esta en el motor.
isFFR = cellfun(@(c) strcmp(c.scheme,'ffr'), CASES);
fprintf('[H2] Delta*SE_centro/SE_borde por ESQUEMA (NaN = no aplica): %s\n', ...
    mat2str(round(OUT.diag.crit_ratio_ffr,3)));
fprintf('     solo cuentan los esquemas FFR (%s): %s -> interior rentable? %s\n', ...
    strjoin(nm(isFFR),', '), mat2str(round(OUT.diag.crit_ratio_ffr(isFFR),3)), ...
    tern(any(OUT.diag.crit_ratio_ffr(isFFR) > 1),'SI','NO'));

%% ============ 4. E3b EQUIVALENTE: BARRIDO DE MASCARA ====================
fprintf('\n=======================================================================\n');
fprintf('  4. EJE DE ELEVACION (equivalente a E3b)\n');
fprintf('=======================================================================\n');
fprintf(['NOTA: 55 deg es el diseno NOMINAL de OneWeb; 45-50 solo se permite en bajas\n' ...
         'latitudes [TN p.13 §A.4]. Los puntos < 55 quedan FUERA del diseno nominal y se\n' ...
         'marcan como tales. Cada mascara reconstruye la tabla P.618 (~80 s).\n\n']);

ELIST = [70 65 60 55 50 45];
RE = cell(1, numel(ELIST));
for i = 1:numel(ELIST)
    pti = struct('cases',{CASES}, 'label',sprintf('minElev=%d',ELIST(i)), ...
                 'cfg_over',struct('geom',struct('minElev',ELIST(i))), 'idx',10+i);
    fprintf('  [%d/%d] minElev = %2d deg ... ', i, numel(ELIST), ELIST(i));
    tt = tic;  RE{i} = run_one_density(cfg, pti);  fprintf('%.0f s\n', toc(tt));
end

fprintf('\n%9s %8s %10s', 'minElev', 'cob.geo', 'elev_med');
for c = 1:numel(nm), fprintf(' %13s', nm{c}); end
fprintf('   (SINR_edge_p5 [dB])\n');
for i = 1:numel(ELIST)
    fprintf('%9d %8.3f %10.1f', ELIST(i), RE{i}.diag.covFrac, RE{i}.diag.elev_mean);
    for c = 1:numel(nm), fprintf(' %13.3f', RE{i}.K{c}.viab.SINR_edge_p5); end
    fprintf('%s\n', tern(ELIST(i) < 55, '   [fuera del nominal]', ''));
end

fprintf('\n%9s', 'minElev');
for c = 1:numel(nm), fprintf(' %13s', nm{c}); end
fprintf('   (veredicto)\n');
for i = 1:numel(ELIST)
    fprintf('%9d', ELIST(i));
    for c = 1:numel(nm), fprintf(' %13s', RE{i}.K{c}.viab.verdict); end
    fprintf('\n');
end

% Margen de viabilidad: mascara MINIMA con verdict = 'viable'
fprintf('\nMARGEN DE VIABILIDAD (mascara MINIMA que sigue siendo ''viable''):\n');
MARG = struct('scheme',nm,'elev_min_viable',num2cell(nan(1,numel(nm))));
for c = 1:numel(nm)
    v = cellfun(@(R) strcmp(R.K{c}.viab.verdict,'viable'), RE);
    if any(v), MARG(c).elev_min_viable = min(ELIST(v)); end
    fprintf('  %-14s %s\n', nm{c}, ...
        tern(any(v), sprintf('%d deg', min(ELIST(v))), 'VACIO (nunca viable)'));
end

%% ============ 5. CONTRASTE DE TOPOLOGIA: 4x4 vs HEXAGONAL ==============
fprintf('\n=======================================================================\n');
fprintf('  5. CONTRASTE DE TOPOLOGIA: reticula 4x4 (16 haces) vs hex (19)\n');
fprintf('=======================================================================\n');
cfgR = config_oneweb('rect');
cfgR.time.dt = cfg.time.dt;  cfgR.time.duration = cfg.time.duration;
cfgR.ffr.tau_mode = 'quantile';  cfgR.ffr.tau_q = TAU_Q;  cfgR.ffr.sched = 'share';
cfgR.compute.timeBlock = cfg.compute.timeBlock;
BLR = build_beam_layout(cfgR);
dR  = sqrt(sum(BLR.cell_km.^2,2));
guardR = (min(dR(BLR.ring == max(BLR.ring))) - cfgR.ground.radius_km)/BLR.spacing_km;
fprintf('layout rect: %d haces | s=%.2f km | extension %.0f km | guarda %.2f anillos\n', ...
    BLR.nBeams, BLR.spacing_km, 2*max(dR), guardR);
fprintf(['OJO: grid_outside_cluster usa el circunradio HEXAGONAL s/sqrt(3)=0.577s. En una\n' ...
         'celda CUADRADA el punto mas lejano es la esquina, a s*sqrt(2)/2=0.707s, luego el\n' ...
         'criterio operativo es CONSERVADOR aqui (avisaria antes de tiempo).\n']);
for D = [3 4]
    [~, COLr] = ffr_coloring(BLR, D, false);
    [~, COLh] = ffr_coloring(BL,  D, false);
    fprintf('  Delta=%d: dmin hex = %.3f*s | dmin rect = %.3f*s  -> %s\n', D, ...
        COLh.dmin_s, COLr.dmin_s, tern(COLr.dmin_s < COLh.dmin_s, 'la cuadrada SEPARA MENOS', 'igual o mejor'));
end
ptR = struct('cases',{CASES}, 'label','OneWeb 4x4', 'idx',30);
tic; OUTR = run_one_density(cfgR, ptR); fprintf('\n[rect] %.0f s | cobertura %.3f\n\n', toc, OUTR.diag.covFrac);
print_scheme_table(OUTR.K, nm, 'TOPOLOGIA 4x4 (16 haces)');

fprintf('\n%-14s %14s %14s %10s   %-10s -> %-10s\n', ...
    'esquema','hex(19)','rect(16)','dif[dB]','ver.hex','ver.rect');
for c = 1:numel(nm)
    a = OUT.K{c}.viab;  b = OUTR.K{c}.viab;
    fprintf('%-14s %14.3f %14.3f %10.3f   %-10s -> %-10s %s\n', nm{c}, ...
        a.SINR_edge_p5, b.SINR_edge_p5, b.SINR_edge_p5-a.SINR_edge_p5, ...
        a.verdict, b.verdict, tern(~strcmp(a.verdict,b.verdict),'*** CAMBIA ***',''));
end

%% ====== 6. E3 EQUIVALENTE: DENSIFICAR HACES CON LA HUELLA FIJA =========
fprintf('\n=======================================================================\n');
fprintf('  6. DENSIFICACION DE HACES CON LA HUELLA FIJA (equivalente a E3)\n');
fprintf('=======================================================================\n');
fprintf(['El eje fisicamente sensato en OneWeb NO es anadir anillos (eso saldria de la\n' ...
         'huella y invadiria celdas de otros satelites, que no son co-canal): es meter MAS\n' ...
         'haces MAS ESTRECHOS dentro de la MISMA huella de 1140 km. Con R anillos y huella\n' ...
         'fija, s = 570/R km. Es lo que haria un OneWeb de nueva generacion.\n\n']);

RLIST = [2 3 4 5];
RD = cell(1,numel(RLIST));  densB = nan(1,numel(RLIST));
fprintf('%8s %8s %10s %11s %16s\n','nRings','nBeams','s[km]','HPBW[deg]','haces/1000km2');
for i = 1:numel(RLIST)
    % SEGUNDO uso INVERSO de la relacion normativa (el primero esta en
    % config_oneweb): se deduce la HPBW de la separacion, con la forma APROXIMADA,
    % la misma que build_beam_layout aplica hacia adelante -> las dos se cancelan y
    % s_km sale exactamente 570/R. Ver el aviso de build_beam_layout: la inversa
    % exacta daria un 0.91% menos, y si alguna vez se cambia hay que cambiar LAS DOS
    % formas a la vez.
    sR   = 570/RLIST(i);
    sdR  = atand(sR/cfg.constellations(1).h);
    hpR  = 2*asind(deg2rad(sdR)/sqrt(3));
    nB   = 3*RLIST(i)*(RLIST(i)+1)+1;
    densB(i) = nB / (1140^2/1000);
    fprintf('%8d %8d %10.1f %11.4f %16.4f\n', RLIST(i), nB, sR, hpR, densB(i));
    ov = struct('radio',struct('beamwidth3dB_deg',hpR), 'beams',struct('nRings',RLIST(i),'nBeams',[]));
    pti = struct('cases',{CASES}, 'label',sprintf('nRings=%d',RLIST(i)), ...
                 'cfg_over',ov, 'idx',40+i);
    RD{i} = run_one_density(cfg, pti);
end
fprintf('\n%8s %11s', 'nRings', 'dens');
for c = 1:numel(nm), fprintf(' %13s', nm{c}); end
fprintf('   (SINR_edge_p5 [dB])\n');
for i = 1:numel(RLIST)
    fprintf('%8d %11.4f', RLIST(i), densB(i));
    for c = 1:numel(nm), fprintf(' %13.3f', RD{i}.K{c}.viab.SINR_edge_p5); end
    fprintf('\n');
end

%% ==================== 7. FIGURAS Y GUARDADO ============================
figure('Name','OneWeb | auto-similitud','Color','w');
semilogx(AS.HPBW, AS.dev, 'o-','LineWidth',1.6); grid on; hold on;
xline(2.08,'--','Starlink'); xline(cfg.radio.beamwidth3dB_deg,'--','OneWeb');
xlabel('HPBW [deg]'); ylabel('desviacion de G(vecino) vs asintotico [dB]');
title('Dominio de validez de la invariancia de escala de la reticula co-canal');
save_fig(gcf, 'figs_oneweb', 'autosimilitud');

figure('Name','OneWeb | eje de elevacion','Color','w');
hold on;
for c = [1 2 4]
    plot(ELIST, cellfun(@(R) R.K{c}.viab.SINR_edge_p5, RE), 'o-','LineWidth',1.6,'DisplayName',nm{c});
end
yline(0,'k--','umbral 0 dB'); yline(-6.7,'r--','suelo QPSK');
xline(55,':','nominal 55 deg');
set(gca,'XDir','reverse'); grid on; legend('Location','best');
xlabel('mascara de elevacion [deg]'); ylabel('SINR_{edge,p5} [dB]');
title('OneWeb: viabilidad de borde frente a la mascara de elevacion');
save_fig(gcf, 'figs_oneweb', 'elevacion');

figure('Name','OneWeb | densificacion','Color','w');
hold on;
for c = [1 2 4]
    plot(densB, cellfun(@(R) R.K{c}.viab.SINR_edge_p5, RD), 'o-','LineWidth',1.6,'DisplayName',nm{c});
end
yline(0,'k--'); yline(-6.7,'r--'); grid on; legend('Location','best');
xlabel('densidad de haces [haces/1000 km^2]'); ylabel('SINR_{edge,p5} [dB]');
title('OneWeb: densificacion de haces con la huella FIJA');
save_fig(gcf, 'figs_oneweb', 'densificacion');

% (d) CONTRASTE DE TOPOLOGIA: hexagonal (19 haces) vs cuadrada 4x4 (16).
%  Es el hallazgo "la ventaja de Delta=3 es un constructo HEXAGONAL" y hasta ahora
%  solo existia como tabla. La clave que explica el grafico es la SEPARACION
%  co-canal: Delta=3 cae de sqrt(3)*s a sqrt(2)*s al pasar a reticula cuadrada,
%  mientras Delta=4 se queda en 2*s en las dos.
figure('Name','OneWeb | topologia hex vs rect','Color','w');
selT = [1 2 3];                                   % reuse1, reuseD(3), reuseD(4)
if numel(nm) >= 5, selT = [1 2 3 4 5]; end        % + ffr(3), ffr(4) si existen
vh = arrayfun(@(c) OUT.K{c}.viab.SINR_edge_p5,  selT);
vr = arrayfun(@(c) OUTR.K{c}.viab.SINR_edge_p5, selT);
hb = bar([vh(:) vr(:)]);
hb(1).FaceColor = [0.20 0.45 0.70];
hb(2).FaceColor = [0.85 0.45 0.15];
hold on;
yline(0,    'k--','LineWidth',1.5,'Label','umbral 0 dB (viable)');
yline(-6.7, 'r--','LineWidth',1.5,'Label','suelo QPSK (-6.7 dB)');
% Etiquetas con la diferencia, que es el resultado
for i = 1:numel(selT)
    d = vr(i) - vh(i);
    text(i, max(vh(i),vr(i)) + 0.6, sprintf('%+.2f dB', d), ...
        'HorizontalAlignment','center', 'FontWeight','bold', ...
        'Color', [0.6 0 0]*(d<0) + [0 0.45 0]*(d>=0));
end
set(gca,'XTickLabel',nm(selT),'XTickLabelRotation',15);
grid on; ylabel('SINR_{edge,p5} [dB]');
legend({'hexagonal (19 haces)','cuadrada 4x4 (16 haces)'}, 'Location','northwest');
title({'La ventaja de \Delta=3 es un constructo HEXAGONAL', ...
       '\Delta=3 pierde (co-canal a \surd2\cdots en vez de \surd3\cdots); \Delta=4 gana (2\cdots en ambas)'});
save_fig(gcf, 'figs_oneweb', 'topologia_hex_vs_rect');

fprintf('\n[figuras] 4 PNG a 300 dpi en figs_oneweb%s\n', filesep);

save('oneweb_results.mat','cfg','cfgR','BL','BLR','CASES','nm','OUT','OUTR', ...
     'RE','ELIST','MARG','RD','RLIST','densB','AS','ALPHA_STAR','TAU_Q');
fprintf('\nResultados en oneweb_results.mat | figuras en figs_oneweb/\n');

%% ========================= FUNCIONES ===================================
function print_scheme_table(K, nm, titulo)
fprintf('\n--- %s ---\n', titulo);
fprintf('%-14s %10s %10s %10s %10s %9s %8s %10s %s\n', ...
    'esquema','SINRe_p5','Rp5_bord','Rp5_glob','R_agg[Gb]','cob.cal','Jain','penaliz','veredicto');
for c = 1:numel(K)
    v = K{c}.viab;
    fprintf('%-14s %10.3f %10.2f %10.2f %10.2f %9.2f %8.3f %10.2f %s\n', nm{c}, ...
        v.SINR_edge_p5, K{c}.edge.R_p5_Mbps, K{c}.all.R_p5_Mbps, ...
        K{c}.all.R_agg_Mbps/1e3, v.coverage, K{c}.all.jain, K{c}.penalty_mean, v.verdict);
end
end

function s = tern(c,a,b)
if c, s = a; else, s = b; end
end
