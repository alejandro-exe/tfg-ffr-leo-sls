function B_Hz = radio_band_Hz(cfg)
%RADIO_BAND_HZ  Ancho de banda por haz en Hz, DERIVADO en el punto de uso.
%   B_Hz = radio_band_Hz(cfg)
%
%   FUENTE UNICA de la conversion MHz -> Hz (misma politica que merge_cfg,
%   mem_peak_model_GB y grid_outside_cluster).
%
%   POR QUE EXISTE (auditoria, hallazgo G3)
%   ---------------------------------------
%   Antes, config_default guardaba el valor YA CALCULADO:
%
%       cfg.radio.B_MHz = 250;
%       cfg.radio.B_Hz  = cfg.radio.B_MHz * 1e6;     % <-- CAMPO DERIVADO, CONGELADO
%
%   El problema no es la formula, es que el resultado quedaba PETRIFICADO en la cfg.
%   Cualquiera que escribiera despues
%
%       cfg.radio.B_MHz = 100;                       % o via pt.cfg_over / merge_cfg
%
%   dejaba B_Hz valiendo 2.5e8 (verificado en la auditoria). Y B_Hz alimenta el RUIDO
%   TERMICO (N = -228.6 + 10log10(Tsys) + 10log10(B_Hz)), luego un barrido de ancho de
%   banda habria salido con el ruido del valor viejo: numeros perfectamente plausibles
%   y silenciosamente equivocados.
%
%   Es EXACTAMENTE el modo de fallo que ya mordio dos veces en este proyecto:
%     - atm_key_of hacia un merge SUPERFICIAL y perdia los campos hermanos de cfg.radio
%     - associate_serving elegia servidor de CUALQUIER constelacion
%   y `pt.cfg_over` existe justamente para barrer hojas arbitrarias de cfg.radio, luego
%   la trampa estaba armada.
%
%   SOLUCION: el campo cfg.radio.B_Hz **ya no existe**. B_MHz es el unico maestro y la
%   conversion se hace aqui, en el punto de uso, cada vez. Asi es ESTRUCTURALMENTE
%   IMPOSIBLE que quede obsoleto: no hay nada que pueda quedarse desactualizado.
%
%   Coste: se llama O(bloques temporales), no O(elementos) -- una vez por
%   compute_link_budget / compute_interference / ffr_allocate. Irrelevante.
%
%   BIT-IDENTICO al comportamiento anterior: 250*1e6 = 2.5e8 es exacto en IEEE-754 y
%   es literalmente la misma operacion que hacia config_default.
%
%   Compatibilidad con cfg antiguas recuperadas de un .mat: si la cfg todavia arrastra
%   un campo B_Hz, assert_cfg_coherent comprueba que coincide y ABORTA si no.
%
%   Aportacion propia (infraestructura; no interviene en ningun calculo fisico).

B_Hz = cfg.radio.B_MHz * 1e6;
end
