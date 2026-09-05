function FSPL_dB = freespace_loss(range_km, freq_GHz)
%FREESPACE_LOSS  Perdida en espacio libre (formula de Friis), programada a mano.
%   FSPL_dB = freespace_loss(range_km, freq_GHz)
%
%       FSPL[dB] = 20*log10(4*pi*d[m]*f[Hz]/c[m/s])
%
%   Aportacion propia (no toolbox). Acepta range_km escalar o matricial; propaga
%   NaN donde el rango es NaN (sin cobertura). freq_GHz escalar.

c    = 299792458;            % m/s
f_Hz = freq_GHz * 1e9;
d_m  = range_km * 1e3;       % km -> m
FSPL_dB = 20 * log10(4*pi .* d_m .* f_Hz / c);
FSPL_dB(isnan(range_km)) = NaN;
end
