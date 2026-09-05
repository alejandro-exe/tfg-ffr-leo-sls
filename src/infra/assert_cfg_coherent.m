function assert_cfg_coherent(cfg)
%ASSERT_CFG_COHERENT  Comprueba los INVARIANTES de cfg.radio y aborta si se rompen.
%   assert_cfg_coherent(cfg)
%
%   FUENTE UNICA de la comprobacion de coherencia de la configuracion (misma politica
%   que merge_cfg, mem_peak_model_GB, grid_outside_cluster y radio_band_Hz).
%
%   POR QUE EXISTE (auditoria, hallazgo G3)
%   ---------------------------------------
%   config_default guardaba DOS campos derivados ya calculados, que quedaban
%   PETRIFICADOS si alguien sobrescribia sus maestros despues:
%
%       cfg.radio.B_Hz   = B_MHz * 1e6
%       cfg.radio.Tsys_K = 10^((Grx_max_dBi - GT_dBK)/10)
%
%   Medido en la auditoria: poner cfg.radio.GT_dBK = 20 dejaba Tsys_K = 251.19 K
%   cuando deberia valer 79.43 K, es decir **5 dB de error de ruido, silencioso**. Y
%   pt.cfg_over / merge_cfg existen justamente para barrer hojas arbitrarias de
%   cfg.radio, luego la trampa estaba armada. Es el mismo modo de fallo que ya mordio
%   en atm_key_of (merge superficial) y en associate_serving (servidor ajeno).
%
%   SE HAN RESUELTO DE DOS MANERAS DISTINTAS, Y LA RAZON HAY QUE PODER DEFENDERLA
%   ----------------------------------------------------------------------------
%   (a) **B_Hz -> ELIMINADO y derivado en el punto de uso** (radio_band_Hz). Aqui hay
%       un unico maestro posible, B_MHz, y la relacion es una conversion de unidades:
%       se puede recalcular siempre, luego lo correcto es que el campo NO EXISTA. Asi
%       es estructuralmente imposible que quede obsoleto.
%
%   (b) **Tsys_K / GT_dBK -> NO se pueden derivar en el punto de uso.** Los dos
%       perfiles del proyecto eligen maestros OPUESTOS, y los dos tienen razon:
%         - config_default (Ku, estudio):  G/T es el MAESTRO (dato de terminal VSAT)
%                                          -> Tsys_K = 10^((Grx_max - GT)/10)
%         - config_calib  (Ka, E0 38.821): Tsys_K es el MAESTRO (Tant + T0*(F-1),
%                                          Tabla 6.1.1.1-3) -> GT_dBK = Grx_max - 10log10(Tsys)
%       Forzar una unica direccion romperia el otro perfil: derivar Tsys de GT en E0
%       obligaria a un ida y vuelta 10^(log10(x)) que perturba el resultado a nivel de
%       precision de maquina (~1e-13) y E0 dejaria de ser bit-identica. Por eso aqui
%       NO se deriva: se COMPRUEBA el invariante, que es el mismo en ambos sentidos:
%
%           GT_dBK == Grx_max_dBi - 10*log10(Tsys_K)
%
%       Si alguien sobrescribe uno de los dos sin recalcular el otro, esto ABORTA en
%       lugar de devolver un ruido termico equivocado con pinta de correcto.
%
%   DONDE SE LLAMA
%     - compute_link_budget : punto de paso UNIVERSAL (todo pipeline calcula balance
%       de enlace: E0, run_ffr_demo, E3, E3b, E4, E5, barridos). Coste O(bloques
%       temporales), no O(elementos): irrelevante.
%     - run_one_density     : ademas al principio de cada punto de barrido, para
%       FALLAR PRONTO -- es donde se aplica pt.cfg_over, y asi un barrido mal
%       configurado aborta antes de gastar minutos de calculo.
%
%   No cambia ningun numero: solo puede abortar donde antes se habria seguido con una
%   cfg incoherente.
%
%   Aportacion propia (infraestructura; no interviene en ningun calculo fisico).

R = cfg.radio;

%% 1. Ancho de banda: maestro unico B_MHz
if ~isfield(R,'B_MHz') || ~isscalar(R.B_MHz) || ~isfinite(R.B_MHz) || R.B_MHz <= 0
    error('assert_cfg_coherent:B_MHz', ...
        'cfg.radio.B_MHz debe ser un escalar finito > 0 (recibido: %s).', mat2str(R.B_MHz));
end

% Compatibilidad: una cfg antigua recuperada de un .mat puede arrastrar todavia el
% campo derivado B_Hz. Si esta y NO coincide, es exactamente el bug que se persigue.
if isfield(R,'B_Hz')
    B_ok = radio_band_Hz(cfg);
    if abs(R.B_Hz - B_ok) > 1e-6 * max(B_ok, 1)
        error('assert_cfg_coherent:B_Hz_stale', ...
            ['cfg.radio.B_Hz = %.6g Hz es INCOHERENTE con cfg.radio.B_MHz = %.6g MHz ' ...
             '(deberia ser %.6g Hz). Es un campo derivado OBSOLETO: B_Hz ya no forma ' ...
             'parte de la configuracion; usa radio_band_Hz(cfg).'], R.B_Hz, R.B_MHz, B_ok);
    end
end

%% 2. Cadena de ruido: G/T, Grx_max y Tsys deben ser mutuamente coherentes
for f = {'Grx_max_dBi','GT_dBK','Tsys_K'}
    if ~isfield(R, f{1}) || ~isscalar(R.(f{1})) || ~isfinite(R.(f{1}))
        error('assert_cfg_coherent:missing', ...
            'cfg.radio.%s ausente o no escalar finito.', f{1});
    end
end
if R.Tsys_K <= 0
    error('assert_cfg_coherent:Tsys', ...
        'cfg.radio.Tsys_K debe ser > 0 K (recibido %.6g).', R.Tsys_K);
end

GT_from_Tsys = R.Grx_max_dBi - 10*log10(R.Tsys_K);
dGT = abs(R.GT_dBK - GT_from_Tsys);
tol = 1e-6;                                  % dB; holgado frente al ida/vuelta 10^log10
if dGT > tol
    error('assert_cfg_coherent:GT_Tsys', ...
        ['INCOHERENCIA en la cadena de ruido (desvio %.4g dB > %.0e):\n' ...
         '  Grx_max_dBi = %.4f dBi\n' ...
         '  GT_dBK      = %.4f dB/K      (declarado)\n' ...
         '  Tsys_K      = %.4f K         -> implica G/T = %.4f dB/K\n' ...
         'Se ha sobrescrito uno sin recalcular el otro. Recuerda cual es el maestro ' ...
         'en cada perfil: config_default -> G/T maestro, Tsys derivado; config_calib ' ...
         '(E0) -> Tsys maestro (Tant + T0*(F-1)), G/T derivado.'], ...
        dGT, tol, R.Grx_max_dBi, R.GT_dBK, R.Tsys_K, GT_from_Tsys);
end
end
