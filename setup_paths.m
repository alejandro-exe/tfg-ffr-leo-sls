function setup_paths()
%SETUP_PATHS  Anade al path de MATLAB las carpetas de codigo del proyecto.
%
%   Ejecutar UNA VEZ al abrir MATLAB, antes de lanzar cualquier experimento:
%
%       setup_paths
%
%   IMPORTANTE: los runners guardan sus resultados con nombres desnudos
%   (save('ffr_results.mat',...)) y sus figuras en carpetas relativas
%   (figs_e2/, figs_e5/...), luego escriben en el DIRECTORIO DE TRABAJO.
%   Ejecuta siempre desde la raiz del proyecto para que todo caiga en el
%   mismo sitio y los anclajes entre experimentos se encuentren.
%
%   Ver el apartado "Puesta en marcha" del README.

HERE = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(HERE, 'src')));
addpath(genpath(fullfile(HERE, 'experimentos')));
addpath(genpath(fullfile(HERE, 'tests')));
addpath(genpath(fullfile(HERE, 'figuras')));

fprintf('Path configurado desde %s\n', HERE);
end
