function gb = mem_peak_model_GB(M, N, Ntb, nB)
%MEM_PEAK_MODEL_GB  Pico de memoria de UN bloque temporal de Ntb instantes [GB].
%   gb = mem_peak_model_GB(M, N, Ntb, nB)
%
%   FUENTE UNICA del modelo de memoria del pipeline (mismo principio que el motor
%   de la FFR: una sola implementacion, reutilizada). La usan run_one_density
%   (para informar OUT.mem del punto que acaba de ejecutar) y estimate_sweep_memory
%   (para dimensionar el barrido ANTES de lanzarlo). Si las dos tuvieran su propia
%   copia, el estimador y lo medido podrian divergir sin que nadie lo notara, que es
%   exactamente el fallo que tenia el estimador de run_sweep_points (25 B/elem, 3x
%   por debajo del modelo real).
%
%   Entradas:
%     M   - usuarios de la rejilla
%     N   - satelites (TODAS las constelaciones sumadas)
%     Ntb - instantes DEL BLOQUE (no de la ventana): con cfg.compute.timeBlock = n
%           el pico lo fija n, no Nt. Sin trocear, Ntb = Nt.
%     nB  - haces del cluster (BL.nBeams)
%
%   MODELO (medido en la FASE A, run_convergence_studyA). El maximo NO son los
%   tensores de compute_geometry sino el interior de compute_interference, en la
%   llamada atm_loss_dB(cfg, G.el) sobre un array [M x N x Ntb], donde coexisten
%   por ELEMENTO:
%       G.el/az/range 24 B + G.vis 1 B + FSPL_all 8 B + Latm_clear/rain 16 B
%       + valid 1 B + temporales de interp1 ~24 B   =   ~74 B/elemento
%   A eso se suma ffr_context, que es M*nBeams*Ntb*(psi + Grel + Glin + temporales)
%   ~40 B/elemento. Contrastado con el pico REAL del proceso: 61 B/elem medidos
%   frente a 74 modelados, luego el modelo es COTA SUPERIOR (que es lo que se
%   quiere para decidir cuantos workers caben).
%
%   Los acumuladores de muestras [M x Nt] del troceado (unos cientos de kB) son
%   despreciables frente a esto y NO entran en el modelo.
%
%   Aportacion propia (infraestructura; no interviene en ningun calculo fisico).

gb = M .* Ntb .* (N .* 74 + nB .* 40) / 2^30;
end
