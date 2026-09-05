%% TEST_ATM_CACHE  Verificacion de la cache (persistent) de atm_loss_dB.
%  La cache es una optimizacion de RENDIMIENTO: NO debe cambiar ningun resultado
%  numerico. Este script lo comprueba en cuatro pruebas:
%
%    (a) resultado de referencia con la cache vacia,
%    (b) INVARIANCIA: tras vaciar la cache el resultado debe ser IDENTICO,
%    (c) RENDIMIENTO: la 2a llamada (tabla cacheada) debe ser mucho mas rapida,
%    (d) VALIDEZ DE LA CLAVE: al cambiar la configuracion la tabla SE RECONSTRUYE
%        (se contrasta contra una llamada DIRECTA a p618PropagationLosses) y al
%        volver a la configuracion original se recupera el resultado exacto.
%
%  Uso: ejecutar en la raiz del proyecto.  No escribe ningun fichero.

clear; clc;
cfg = config_default();
ok  = true;

fprintf('===== TEST DE LA CACHE DE atm_loss_dB =====\n');
fprintf('cfg: f=%.4g GHz | punto=(%.4f, %.4f) | minElev=%g deg | P.618 exc=%g%%\n\n', ...
    cfg.radio.freq_GHz, cfg.ground.point(1), cfg.ground.point(2), ...
    cfg.geom.minElev, cfg.radio.p618_availability);

% Elevaciones variadas: nodos exactos de la tabla, valores intermedios, el limite
% superior, un valor POR DEBAJO de minElev (extrapolacion) y NaN (sin cobertura).
el = [cfg.geom.minElev, 25.7, 27, 33.3, 35, 47.5, 60, 72.4, 89, 89.9, 5, NaN];
elMat = reshape([30 40 50 60 70 80], 2, 3);      % comprueba que respeta el tamano

%% (a) Referencia con la cache VACIA
clear atm_loss_dB
[c1, r1, u1] = atm_loss_dB(cfg, el);
[cm1, rm1]   = atm_loss_dB(cfg, elMat);
fprintf('(a) referencia: p618_used=%d | Latm_clear(35 deg)=%.6f dB | Latm_rain(35 deg)=%.6f dB\n', ...
    u1, c1(5), r1(5));

%% (b) INVARIANCIA: vaciar la cache y repetir -> identico bit a bit
clear atm_loss_dB
[c2, r2, u2] = atm_loss_dB(cfg, el);
[cm2, rm2]   = atm_loss_dB(cfg, elMat);

% El vector incluye NaN A PROPOSITO (elevacion sin cobertura), e isequal devuelve
% false ante NaN==NaN. Por eso se comprueban las dos cosas: isequal EXACTO sobre
% las muestras no-NaN, e isequaln (que trata NaN como igual) sobre el array entero.
fin  = ~isnan(el);
eqEx = isequal(c1(fin), c2(fin)) && isequal(r1(fin), r2(fin));
eqAll= isequaln(c1, c2) && isequaln(r1, r2) && isequal(u1, u2) && ...
       isequaln(cm1, cm2) && isequaln(rm1, rm2);
fprintf('(b) invariancia tras vaciar la cache: isequal (no-NaN) = %d | isequaln (todo) = %d | max|dif| = %g\n', ...
    eqEx, eqAll, max([abs(c1(fin)-c2(fin)), abs(r1(fin)-r2(fin))]));
if ~(eqEx && eqAll)
    ok = false; fprintf('    ==> FALLO: la cache ALTERA el resultado.\n');
else
    fprintf('    ==> OK: resultado IDENTICO bit a bit.\n');
end
if ~isequal(size(cm1), size(elMat)), ok = false; fprintf('    ==> FALLO: no respeta el tamano de entrada.\n'); end

%% (c) RENDIMIENTO: 1a llamada (construye tabla) vs 2a (tabla cacheada)
elBig = linspace(cfg.geom.minElev, 89, 5000);
clear atm_loss_dB
t0 = tic;  atm_loss_dB(cfg, elBig);  t_cold = toc(t0);   % construye la tabla
t1 = tic;  atm_loss_dB(cfg, elBig);  t_warm = toc(t1);   % ya cacheada

nRep = 20;  tw = zeros(1, nRep);
for i = 1:nRep, t = tic; atm_loss_dB(cfg, elBig); tw(i) = toc(t); end
t_warm_med = median(tw);

fprintf('\n(c) rendimiento (%d elevaciones):\n', numel(elBig));
fprintf('    1a llamada (construye la tabla) : %8.2f ms\n', 1e3*t_cold);
fprintf('    2a llamada (tabla cacheada)     : %8.3f ms   -> speedup x%.0f\n', 1e3*t_warm, t_cold/t_warm);
fprintf('    mediana de %2d llamadas calientes: %8.3f ms   -> speedup x%.0f\n', ...
    nRep, 1e3*t_warm_med, t_cold/t_warm_med);
if t_warm_med >= t_cold
    ok = false; fprintf('    ==> FALLO: la 2a llamada no es mas rapida (la cache no actua).\n');
else
    fprintf('    ==> OK: la tabla se construye UNA sola vez.\n');
end

%% (d) VALIDEZ DE LA CLAVE: cambiar freq_GHz debe RECONSTRUIR la tabla
elNode = cfg.geom.minElev + 10;      % nodo EXACTO de la tabla (paso 2 deg)
cfg2 = cfg;  cfg2.radio.freq_GHz = 20;      % Ku -> Ka
[c3, r3] = atm_loss_dB(cfg2, elNode);

% Calculo DIRECTO conocido (lo que deberia devolver a 20 GHz), sin pasar por la cache
cCS = p618Config('Frequency', cfg2.radio.freq_GHz*1e9, 'ElevationAngle', elNode, ...
    'Latitude', cfg2.ground.point(1), 'Longitude', cfg2.ground.point(2), ...
    'TotalAnnualExceedance', 50);
cRN = p618Config('Frequency', cfg2.radio.freq_GHz*1e9, 'ElevationAngle', elNode, ...
    'Latitude', cfg2.ground.point(1), 'Longitude', cfg2.ground.point(2), ...
    'TotalAnnualExceedance', cfg2.radio.p618_availability);
resCS = p618PropagationLosses(cCS);   cRef = max(0, resCS.At);
resRN = p618PropagationLosses(cRN);   rRef = max(0, resRN.At);

fprintf('\n(d) cambio de configuracion (freq 12 -> 20 GHz) a elev = %g deg:\n', elNode);
fprintf('    cacheada a 12 GHz : Latm_clear = %.6f dB | Latm_rain = %.6f dB\n', c1(5), r1(5));
fprintf('    devuelto a 20 GHz : Latm_clear = %.6f dB | Latm_rain = %.6f dB\n', c3, r3);
fprintf('    directo  a 20 GHz : Latm_clear = %.6f dB | Latm_rain = %.6f dB\n', cRef, rRef);
if ~(isequal(c3, cRef) && isequal(r3, rRef))
    ok = false; fprintf('    ==> FALLO: no coincide con el calculo directo (tabla obsoleta o mal reconstruida).\n');
elseif isequal(c3, c1(5))
    ok = false; fprintf('    ==> FALLO: devuelve el valor de 12 GHz (la clave no detecta el cambio).\n');
else
    fprintf('    ==> OK: la tabla SE RECONSTRUYE y coincide EXACTAMENTE con el calculo directo.\n');
end

% Round-trip: volver a la configuracion original debe reproducir (a) exactamente
[c4, r4, u4] = atm_loss_dB(cfg, el);
rt = isequaln(c4, c1) && isequaln(r4, r1) && isequal(u4, u1);
fprintf('    round-trip 20 -> 12 GHz reproduce la referencia: %d\n', rt);
if ~rt, ok = false; fprintf('    ==> FALLO: el round-trip no recupera el resultado original.\n'); end

% Los OTROS campos de la clave tambien deben disparar la reconstruccion
pertName = {'radio.freq_GHz', 'ground.point(1) [lat]', 'ground.point(2) [lon]', ...
            'geom.minElev',   'radio.p618_availability'};
fprintf('\n    campos de la clave que deben invalidar la cache:\n');
for i = 1:numel(pertName)
    cfgP = cfg;
    switch i
        case 1, cfgP.radio.freq_GHz          = cfg.radio.freq_GHz + 8;
        case 2, cfgP.ground.point(1)         = cfg.ground.point(1) + 5;
        case 3, cfgP.ground.point(2)         = cfg.ground.point(2) + 5;
        case 4, cfgP.geom.minElev            = cfg.geom.minElev + 1;
        case 5, cfgP.radio.p618_availability = 1.0;
    end
    [cp, rp] = atm_loss_dB(cfgP, elNode);
    [cb, rb] = atm_loss_dB(cfg,  elNode);          % vuelve a la original
    % Nota: p618_availability solo afecta a la rama de LLUVIA (el cielo claro usa
    % siempre exceedance 50%), luego para ese campo el cambio se ve en Latm_rain.
    changed = ~isequal(cp, cb) || ~isequal(rp, rb);
    fprintf('      %-26s -> reconstruye: %d\n', pertName{i}, changed);
    if ~changed, ok = false; fprintf('        ==> FALLO: el cambio NO invalida la cache.\n'); end
end

%% Veredicto
fprintf('\n=====================================================\n');
if ok
    fprintf('TODAS LAS COMPROBACIONES PASAN: la cache no altera ningun resultado.\n');
else
    error('test_atm_cache:fail', 'ALGUNA COMPROBACION HA FALLADO (ver arriba).');
end
fprintf('=====================================================\n');
