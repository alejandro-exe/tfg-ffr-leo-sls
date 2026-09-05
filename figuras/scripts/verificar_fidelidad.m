%% VERIFICAR_FIDELIDAD
%  Comprueba que las figuras REGENERADAS reproducen las publicadas.
%  Solo LEE los PNG publicados; escribe unicamente dentro de REDACCION/.
%
%  POR QUE SE COMPARA LA "REPLICA" Y NO LA VERSION SIN TITULO.
%  save_fig usa exportgraphics con recorte ajustado al CONTENIDO, asi que quitar el
%  title() no se limita a borrar una banda: cambia el lienzo (en fig 4.5 el titulo
%  es MAS ANCHO que los ejes, luego al quitarlo la imagen sale mas estrecha). Una
%  comparacion pixel a pixel sin titulo seria imposible de interpretar. Por eso
%  cada script de regeneracion emite tambien una REPLICA con el titulo original:
%  si la replica coincide con el PNG publicado, queda demostrado que el .mat SI
%  corresponde a la figura publicada, y la version sin titulo difiere de ella
%  UNICAMENTE en la llamada a title() que se ha omitido.
%
%  Uso: desde cualquier sitio,  run('.../figuras/scripts/verificar_fidelidad.m')

HERE     = fileparts(mfilename('fullpath'));
PROJROOT = fileparts(fileparts(HERE));
VERDIR   = fullfile(PROJROOT,'REDACCION','figuras','_verificacion');

pares = {
  fullfile(PROJROOT,'figs_geom','geom_f_constelacion_3d_ecef_t0.png'), ...
  fullfile(VERDIR,'fig_4_02_replica.png'), 'fig 4.2 (constelacion ECEF)'
  fullfile(PROJROOT,'figs_e2','e2_e_celdas_earthfixed_coloreado.png'), ...
  fullfile(VERDIR,'fig_4_05_replica.png'), 'fig 4.5 (celdas coloreadas)'
};

fprintf('\n===== FIDELIDAD: PNG publicado vs REPLICA regenerada (con titulo) =====\n');

% Renderizador en uso: es la explicacion candidata de cualquier diferencia de
% antialiasing. `matlab -batch` corre sin display y suele caer en OpenGL SOFTWARE,
% mientras que la figura publicada se genero en sesion interactiva (hardware).
try
    gl = opengl('data');
    fprintf('\n[renderer] Version : %s\n', gl.Version);
    fprintf('[renderer] Vendor  : %s\n', gl.Vendor);
    fprintf('[renderer] Software: %d  (1 = OpenGL por software)\n', gl.Software);
catch ME
    fprintf('\n[renderer] no disponible: %s\n', ME.message);
end

for k = 1:size(pares,1)
    fA = pares{k,1};  fB = pares{k,2};  nm = pares{k,3};
    A = imread(fA);   B = imread(fB);

    fprintf('\n--- %s ---\n', nm);
    fprintf('  publicado : %d x %d px\n', size(A,2), size(A,1));
    fprintf('  replica   : %d x %d px\n', size(B,2), size(B,1));

    if ~isequal(size(A), size(B))
        fprintf('  [PARAR] DIMENSIONES DISTINTAS: no es reproducible tal cual. REPORTAR.\n');
        continue
    end

    D    = abs(double(A) - double(B));
    maxd = max(D(:));
    dif  = any(D > 0, 3);
    pdif = 100 * mean(dif, 'all');

    fprintf('  max|dif| por canal : %g (de 255)\n', maxd);
    fprintf('  pixeles distintos  : %.4f %% (%d de %d)\n', pdif, nnz(dif), numel(dif));

    if maxd == 0
        fprintf('  VEREDICTO: BIT-IDENTICAS. Fidelidad CONFIRMADA.\n');
    else
        % Localizar las diferencias para poder diagnosticar (renderizado vs datos)
        [rr, cc] = find(dif);
        fprintf('  bbox de las diferencias : filas %d-%d de %d, cols %d-%d de %d\n', ...
                min(rr), max(rr), size(dif,1), min(cc), max(cc), size(dif,2));
        fprintf('  dif media donde difiere : %.2f niveles\n', mean(D(repmat(dif,1,1,size(D,3)))));

        % DIAGNOSTICO: donde se concentran. Si son lineas de rejilla / borde de caja
        % / un glifo de tick, la diferencia es de RASTERIZADO, no de datos.
        perCol = sum(dif, 1);
        [sc, ic] = sort(perCol, 'descend');
        fprintf('  columnas con mas dif    : %s\n', ...
            strjoin(arrayfun(@(a,b) sprintf('col %d (%d px)', a, b), ...
                    ic(1:3), sc(1:3), 'UniformOutput', false), ', '));
        fprintf('  concentracion           : las 3 peores columnas acumulan %.1f %% de las dif.\n', ...
                100*sum(sc(1:3))/nnz(dif));
        % Cuantas diferencias sobreviven si se ignoran las 5 peores columnas
        keep = true(1, size(dif,2)); keep(ic(1:5)) = false;
        fprintf('  dif fuera de esas 5 col : %d px (%.4f %% del total de la imagen)\n', ...
                nnz(dif(:, keep)), 100*nnz(dif(:, keep))/numel(dif));
        fprintf('  VEREDICTO: hay diferencias. Ver el mapa de diferencias y REPORTAR.\n');

        % Mapa de diferencias (rojo donde difiere) sobre la imagen publicada en gris
        gr  = repmat(uint8(0.6*double(rgb2gray(A)) + 100), 1, 1, 3);
        vis = gr;
        Rch = vis(:,:,1); Gch = vis(:,:,2); Bch = vis(:,:,3);
        Rch(dif) = 255; Gch(dif) = 0; Bch(dif) = 0;
        vis(:,:,1) = Rch; vis(:,:,2) = Gch; vis(:,:,3) = Bch;
        outp = fullfile(VERDIR, sprintf('diffmap_%s.png', regexprep(nm,'[^A-Za-z0-9]','_')));
        imwrite(vis, outp);
        fprintf('  mapa de diferencias     : %s\n', outp);
    end
end

fprintf('\n=======================================================================\n');
