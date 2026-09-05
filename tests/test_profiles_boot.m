%% TEST_PROFILES_BOOT  Los TRES perfiles de cfg arrancan y completan un punto minimo.
%  assert_cfg_coherent verifica los invariantes de cfg.radio, pero eso no basta:
%  hace falta recorrer config_default / config_calib / config_oneweb y comprobar
%  que los tres siguen construyendose y atravesando el pipeline. El modo de fallo
%  que cubre es que un perfil quede roto por un cambio transversal y no se entere
%  nadie hasta que se necesita, porque su .mat ya estaba publicado.
%
%  QUE COMPRUEBA, por perfil:
%    1. la funcion de perfil devuelve una cfg sin error;
%    2. assert_cfg_coherent(cfg) pasa (invariante G/T = Grx_max - 10log10(Tsys),
%       ausencia de B_Hz petrificado, ...);
%    3. un punto MINIMO del pipeline completo se ejecuta de extremo a extremo
%       (build_constellation -> ... -> compute_kpis) y devuelve KPIs finitos y un
%       veredicto de viabilidad valido.
%
%  QUE **NO** COMPRUEBA: los NUMEROS. No es un test de regresion -- para eso estan
%  los .mat publicados y la comparacion bit a bit de E0. Es un test de ARRANQUE:
%  responde "el perfil sigue vivo?", no "el perfil sigue dando lo mismo?".
%
%  COSTE (~4 min): NO es el pipeline, que aqui es trivial (rejilla de 5 usuarios y
%  5 instantes, M*N*Nt ~ 4e4). Es la tabla ITU-R P.618: los tres perfiles tienen
%  CLAVES DE CACHE DISTINTAS (Ku 12 GHz/25 deg, Ka 20 GHz/30 deg, OneWeb 12 GHz/55 deg),
%  luego se construyen TRES tablas de ~80 s. Es irreducible sin falsear el test.
%
%  LA CONSTELACION NO SE REDUCE: cada perfil usa la suya (1584 / su propia / 720).
%  Lo que se encoge es la rejilla y la ventana, que es lo que domina el coste
%  M*N*Nt. Asi el test ejercita la geometria REAL de cada perfil.

clear; clc;

fprintf('==========================================================================\n');
fprintf('  TEST DE ARRANQUE DE LOS TRES PERFILES DE cfg\n');
fprintf('==========================================================================\n');

PROFILES = { 'config_default', @config_default,  'Ku Starlink (estudio)'      ; ...
             'config_calib',   @config_calib,    'Ka Set-1 LEO-600 (E0)'      ; ...
             'config_oneweb',  @config_oneweb,   'Ku OneWeb (validacion cruzada)' };
nP = size(PROFILES,1);

RES  = repmat(struct('name','','built',false,'coherent',false,'pipeline',false, ...
                     'verdict','','SINR_edge_p5',NaN,'covFrac',NaN,'secs',NaN, ...
                     'err',''), 1, nP);
tAll = tic;

for p = 1:nP
    name = PROFILES{p,1};  fh = PROFILES{p,2};  desc = PROFILES{p,3};
    RES(p).name = name;
    fprintf('\n--------------------------------------------------------------------\n');
    fprintf('[%d/%d] %s  --  %s\n', p, nP, name, desc);
    tP = tic;

    try
        % --- 1. Construccion del perfil ------------------------------------
        cfg = fh();
        RES(p).built = true;
        fprintf('  (1) cfg construida ......... OK  (f=%.1f GHz, h=%.0f km, minElev=%d deg, HPBW=%.4f deg)\n', ...
            cfg.radio.freq_GHz, cfg.constellations(1).h, cfg.geom.minElev, ...
            cfg.radio.beamwidth3dB_deg);

        % --- 2. Coherencia de cfg -----------------------------------------
        %  Aborta si se rompe el invariante G/T <-> Tsys o si arrastra un B_Hz
        %  obsoleto. Los dos perfiles eligen MAESTROS OPUESTOS (Ku: G/T maestro;
        %  Ka: Tsys maestro) y el invariante es lo unico comun a los dos.
        assert_cfg_coherent(cfg);
        RES(p).coherent = true;
        fprintf('  (2) assert_cfg_coherent .... OK  (G/T=%.2f dB/K, Tsys=%.1f K, Grx=%.1f dBi)\n', ...
            cfg.radio.GT_dBK, cfg.radio.Tsys_K, cfg.radio.Grx_max_dBi);

        % --- 3. Punto MINIMO del pipeline completo -------------------------
        %  Rejilla y ventana encogidas; constelacion y layout INTACTOS.
        cfg.ground.mode      = 'localgrid';
        cfg.ground.radius_km = 20;
        cfg.ground.step_km   = 20;      % 5 usuarios (centro + 4 vecinos)
        cfg.time.dt          = 600;     % s
        cfg.time.duration    = 2400;    % s -> Nt = 5 instantes
        cfg.ffr.tau_mode     = 'quantile';
        cfg.ffr.tau_q        = 50;

        pt = struct('label', name, ...
                    'cases', {{ struct('scheme','reuse1','Delta',1,'alpha',NaN, ...
                                       'adaptive',false,'name','reuse1'), ...
                                struct('scheme','ffr','Delta',3,'alpha',0.5, ...
                                       'adaptive',false,'name','ffr(D=3)'), ...
                                struct('scheme','ffr','Delta',3,'alpha',0.5, ...
                                       'adaptive',true, 'name','ffr-adaptativa') }}, ...
                    'keep_cdf', false, 'seed', 1);

        OUT = run_one_density(cfg, pt);
        RES(p).pipeline = true;

        V = OUT.K{2}.viab;                       % veredicto del esquema ffr(D=3)
        RES(p).verdict      = V.verdict;
        RES(p).SINR_edge_p5 = V.SINR_edge_p5;
        RES(p).covFrac      = OUT.diag.covFrac;

        % El veredicto tiene que ser uno de los tres niveles definidos: si el
        % motor devolviera otra cosa, el criterio de viabilidad estaria roto.
        if ~ismember(lower(V.verdict), {'inviable','marginal','viable'})
            error('test_profiles_boot:verdict', ...
                'Veredicto no reconocido en %s: "%s"', name, V.verdict);
        end
        fprintf('  (3) pipeline completo ...... OK  (%d usuarios x %d instantes x %d sat, %d haces)\n', ...
            OUT.M, OUT.Nt, OUT.N, OUT.nBeams);
        fprintf('      cobertura %.3f | ffr(D=3): SINR_edge_p5 = %+.2f dB -> %s\n', ...
            OUT.diag.covFrac, V.SINR_edge_p5, upper(V.verdict));

    catch ME
        RES(p).err = sprintf('%s: %s', ME.identifier, ME.message);
        fprintf('  *** FALLO: %s\n', RES(p).err);
    end

    RES(p).secs = toc(tP);
    fprintf('  tiempo: %.1f s\n', RES(p).secs);
end

%% Resumen
fprintf('\n==========================================================================\n');
fprintf('%-16s %8s %10s %10s %12s %10s\n', ...
    'perfil','cfg','coherente','pipeline','veredicto','tiempo');
fprintf('%s\n', repmat('-',1,70));
for p = 1:nP
    fprintf('%-16s %8s %10s %10s %12s %9.1fs\n', RES(p).name, ...
        yn(RES(p).built), yn(RES(p).coherent), yn(RES(p).pipeline), ...
        tern(isempty(RES(p).verdict),'-',upper(RES(p).verdict)), RES(p).secs);
end
fprintf('%s\n', repmat('-',1,70));

okAll = all([RES.built] & [RES.coherent] & [RES.pipeline]);
fprintf('\nTiempo total: %.1f s (%.1f min; ~3 tablas P.618 distintas)\n', ...
    toc(tAll), toc(tAll)/60);
if okAll
    fprintf('RESULTADO: los %d perfiles ARRANCAN y completan el pipeline. TEST SUPERADO.\n', nP);
else
    fprintf('RESULTADO: *** ALGUN PERFIL FALLA *** (ver el detalle arriba)\n');
    error('test_profiles_boot:failed','Al menos un perfil no arranca.');
end

% =========================================================================
function s = yn(b),        if b, s = 'SI'; else, s = 'NO'; end,  end
function s = tern(c,a,b),  if c, s = a;    else, s = b;    end,  end
