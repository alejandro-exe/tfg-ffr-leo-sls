function fpath = save_fig(fig, folder, name, dpi)
%SAVE_FIG  Exporta una figura a PNG con resolucion de memoria. FUENTE UNICA.
%   fpath = SAVE_FIG(fig, folder, name)        usa 300 dpi
%   fpath = SAVE_FIG(fig, folder, name, dpi)   fija la resolucion
%
%   fig    : handle de figura (o [] / gcf para la actual)
%   folder : carpeta de destino (se crea si no existe), p.ej. 'figs_e2'
%   name   : nombre del fichero SIN extension. DESCRIPTIVO, no 'fig1'/'fig2':
%            debe entenderse que muestra la figura sin abrirla.
%   dpi    : resolucion en puntos por pulgada (defecto 300, calidad de imprenta)
%
%   POR QUE UN HELPER Y NO exportgraphics DIRECTO EN CADA RUNNER:
%   es la misma politica de FUENTE UNICA que merge_cfg / grid_outside_cluster /
%   mem_peak_model_GB. La resolucion, el recorte de margenes y el color de fondo
%   quedan fijados en un solo sitio, de modo que las ~25 figuras de la memoria
%   salen homogeneas y cambiar el criterio (p.ej. pasar a 600 dpi o a vectorial)
%   es un cambio de UNA linea y no de 25.
%
%   NOTA: 'ContentType','image' fuerza rasterizado. Es deliberado: varias figuras
%   llevan superficies/parches con transparencia (mapa de celdas, skyplot,
%   constelacion 3D) que en vectorial salen mal en muchos visores de PDF.
%
%   NO afecta a ningun resultado: solo escribe ficheros.

if nargin < 4 || isempty(dpi), dpi = 300; end                      % <-- DECISION
if isempty(fig), fig = gcf; end

if ~exist(folder, 'dir')
    mkdir(folder);
end

fpath = fullfile(folder, [name '.png']);

% Fondo blanco explicito: si la figura se creo sin 'Color','w' el PNG saldria
% con el gris por defecto de MATLAB, que en la memoria impresa se nota.
set(fig, 'Color', 'w', 'InvertHardcopy', 'off');

exportgraphics(fig, fpath, 'Resolution', dpi, 'ContentType', 'image', ...
    'BackgroundColor', 'white');
end
