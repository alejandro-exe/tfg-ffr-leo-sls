function [Latm_clear, Latm_rain, p618_used] = atm_loss_dB(cfg, el_deg)
%ATM_LOSS_DB  Perdidas atmosfericas ITU-R P.618 por interpolacion en elevacion.
%   [Latm_clear, Latm_rain, p618_used] = atm_loss_dB(cfg, el_deg)
%
%   Construye una tabla de perdidas vs elevacion (minElev:2:89 deg) llamando a
%   p618PropagationLosses (Satellite Communications Toolbox) y la interpola en
%   los angulos de elevacion pedidos. Precalcular la tabla es mas eficiente que
%   llamar a la toolbox en cada (usuario, satelite, instante): la funcion es
%   determinista para una ubicacion fija.
%
%   CACHE (persistent) -- optimizacion de RENDIMIENTO, sin efecto numerico:
%   la tabla depende SOLO de la configuracion, no de las elevaciones pedidas, y
%   antes se reconstruia (2*n_el llamadas a la toolbox) en CADA invocacion. En un
%   barrido de saturacion tipo E3 (densidad x Monte Carlo) eso son miles de
%   reconstrucciones identicas. Ahora se construye UNA VEZ por configuracion y se
%   reutiliza mientras no cambie ninguno de los parametros que la determinan:
%
%       freq_GHz | ground.point(1) | ground.point(2) | geom.minElev | p618_availability
%
%   Si cualquiera cambia (p.ej. al pasar del perfil Ku de config_default al perfil
%   Ka de config_calib) la clave no coincide y la tabla SE RECONSTRUYE, de modo que
%   nunca se devuelven valores de otra configuracion. La interpolacion se hace
%   SIEMPRE con la tabla vigente, luego los resultados son bit-identicos a los de
%   la version sin cache (verificado en test_atm_cache.m).
%   Para forzar el vaciado de la cache:  clear atm_loss_dB
%
%   "Cielo claro" = excedencia 50% (lluvia ~0, solo gas + centelleo).
%   "Lluvia"      = excedencia cfg.radio.p618_availability% (p.ej. 0.1% = 99.9%).
%
%   Si los mapas digitales ITU no estan instalados, cae a un modelo simplificado
%   de respaldo (gas 0.05/sin(el), lluvia 3/sin(el)) con aviso. La logica de
%   respaldo es identica a la anterior; el unico cambio de comportamiento (no
%   numerico) es que el warning se emite una vez por RECONSTRUCCION de la tabla en
%   lugar de una vez por llamada.
%
%   Toolbox (no aportacion). Modulo reutilizado por compute_link_budget (servidor)
%   y compute_interference (interferentes).
%
%   Entradas:
%     cfg     - configuracion (usa cfg.geom.minElev, cfg.radio.freq_GHz,
%               cfg.ground.point, cfg.radio.p618_availability)
%     el_deg  - elevaciones [cualquier tamano], NaN permitido (-> NaN en salida)
%   Salidas (mismo tamano que el_deg):
%     Latm_clear, Latm_rain - perdidas atmosfericas [dB]
%     p618_used             - true si la toolbox P.618 estuvo disponible

persistent cacheKey el_tbl Lclear_tbl Lrain_tbl p618_flag

% Clave de configuracion: TODO lo que determina la tabla. Si algo de esto cambia,
% la tabla cacheada es invalida y hay que reconstruirla.
key = sprintf('%.10g|%.10g|%.10g|%.10g|%.10g', ...
    cfg.radio.freq_GHz, cfg.ground.point(1), cfg.ground.point(2), ...
    cfg.geom.minElev,   cfg.radio.p618_availability);

if isempty(cacheKey) || ~strcmp(cacheKey, key)
    %% --- Reconstruccion de la tabla (solo al cambiar la configuracion) ---
    f_Hz    = cfg.radio.freq_GHz * 1e9;
    el_new  = (cfg.geom.minElev : 2 : 89);     % deg (89 evita singularidad a 90)
    n_el    = numel(el_new);
    Lc_new  = zeros(1, n_el);
    Lr_new  = zeros(1, n_el);
    used    = false;

    for i = 1:n_el
        try
            cfgCS = p618Config( ...
                'Frequency',             f_Hz, ...
                'ElevationAngle',        el_new(i), ...
                'Latitude',              cfg.ground.point(1), ...
                'Longitude',             cfg.ground.point(2), ...
                'TotalAnnualExceedance', 50);
            rCS = p618PropagationLosses(cfgCS);
            Lc_new(i) = rCS.At;

            cfgRN = p618Config( ...
                'Frequency',             f_Hz, ...
                'ElevationAngle',        el_new(i), ...
                'Latitude',              cfg.ground.point(1), ...
                'Longitude',             cfg.ground.point(2), ...
                'TotalAnnualExceedance', cfg.radio.p618_availability);
            rRN = p618PropagationLosses(cfgRN);
            Lr_new(i) = rRN.At;

            used = true;
        catch ME
            % Modelo aproximado de respaldo (no toolbox).
            Lc_new(i) = 0.05 / sind(el_new(i));
            Lr_new(i) = 3.0  / sind(el_new(i));
            if i == 1
                warning('atm_loss_dB:p618fallback', ...
                    'p618PropagationLosses no disponible (%s).\nUsando modelo simplificado de respaldo.', ...
                    ME.message);
            end
        end
    end

    % La clave se fija LA ULTIMA: si la construccion se interrumpe (error o Ctrl-C)
    % no queda una clave valida apuntando a una tabla a medio construir.
    el_tbl     = el_new;
    Lclear_tbl = Lc_new;
    Lrain_tbl  = Lr_new;
    p618_flag  = used;
    cacheKey   = key;
end

p618_used = p618_flag;

%% --- Interpolacion: SIEMPRE con la tabla vigente ---
Latm_clear = nan(size(el_deg));
Latm_rain  = nan(size(el_deg));
valid = ~isnan(el_deg);
if any(valid(:))
    Latm_clear(valid) = max(0, interp1(el_tbl, Lclear_tbl, el_deg(valid), 'linear', 'extrap'));
    Latm_rain(valid)  = max(0, interp1(el_tbl, Lrain_tbl,  el_deg(valid), 'linear', 'extrap'));
end
end
