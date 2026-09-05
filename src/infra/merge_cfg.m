function base = merge_cfg(base, over)
%MERGE_CFG  Funde `over` sobre `base` recursivamente (solo pisa las hojas dadas).
%   cfg = merge_cfg(cfg, over)
%
%   FUENTE UNICA de la regla de override de configuracion. La usan run_one_density
%   (para aplicar pt.cfg_over al punto) y run_sweep_points/atm_key_of (para calcular
%   la clave de la cache P.618 del punto). Tener DOS reglas distintas es justo el
%   fallo que hubo aqui: atm_key_of hacia un merge SUPERFICIAL
%   (cfg.radio = over.radio), que REEMPLAZA el struct entero y borra los campos no
%   mencionados -- con pt.cfg_over.radio.beamwidth3dB_deg desaparecian freq_GHz,
%   B_MHz, etc. y la clave de la cache fallaba. Con la regla compartida, lo que se
%   simula y lo que se diagnostica no pueden divergir.
%
%   REGLA: un campo que es struct ESCALAR en AMBOS se recorre campo a campo; todo lo
%   demas (numero, cadena, cell, [], array de structs como cfg.constellations)
%   SUSTITUYE tal cual. Asi
%       over.radio.beamwidth3dB_deg = 0.7
%   cambia SOLO ese campo y deja intacto el resto de cfg.radio, mientras que
%       over.constellations = <array 1x2>
%   reemplaza la lista entera de constelaciones, que es lo que se quiere.
%
%   Nota: `[]` es un valor VALIDO y sustituye (cfg.beams.nBeams = [] es la forma
%   documentada de decirle a build_beam_layout que mande nRings), asi que NO se
%   filtra por isempty.
%
%   Aportacion propia (infraestructura; no interviene en ningun calculo fisico).

if ~isstruct(over) || isempty(fieldnames(over))
    return;
end

f = fieldnames(over);
for i = 1:numel(f)
    k = f{i};
    if isfield(base, k) && isstruct(base.(k)) && isstruct(over.(k)) ...
            && isscalar(base.(k)) && isscalar(over.(k))
        base.(k) = merge_cfg(base.(k), over.(k));
    else
        base.(k) = over.(k);
    end
end
end
