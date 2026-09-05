function S = associate_serving(cfg, G, satMask)
%ASSOCIATE_SERVING  Selecciona el satelite servidor de cada usuario e instante.
%   S = associate_serving(cfg, G)
%   S = associate_serving(cfg, G, satMask)
%
%   Criterio cfg.geom.serving:
%     'maxElev'  -> el satelite visible de mayor elevacion (mejor enlace)
%     'minRange' -> el satelite visible mas cercano
%
%   satMask (OPCIONAL, [N x 1] logico): satelites ELEGIBLES como servidores.
%   Por defecto todos, luego los llamadores existentes NO cambian.
%   IMPRESCINDIBLE en escenarios MULTI-OPERADOR (E4): sin el, el usuario se
%   asociaria al mejor satelite de CUALQUIER constelacion, es decir se conectaria
%   al satelite de otra empresa. Eso no solo es irreal (es cliente de UN operador)
%   sino que rompe el layout Earth-fixed, que se proyecta con la altitud de
%   cfg.constellations(1): un servidor a 1200 km apuntando celdas dimensionadas
%   para 550 km hunde la SINR y el efecto se confundiria con interferencia
%   inter-constelacion. Los satelites ajenos siguen contando como INTERFERENTES
%   (S.nVis los cuenta y compute_interference los usa); solo se les excluye de ser
%   servidores.
%
%   Devuelve:
%     S.idx   [M x Nt]  indice del satelite servidor (NaN si no hay cobertura)
%     S.el    [M x Nt]  elevacion del servidor [deg]
%     S.range [M x Nt]  rango del servidor [km]
%     S.nVis  [M x Nt]  numero de satelites visibles (futuros interferentes; se
%                       cuentan TODOS, elegibles o no)

[M,N,Nt] = size(G.el);
if nargin < 3 || isempty(satMask)
    satMask = true(N,1);
else
    satMask = logical(satMask(:));
end
S.idx   = nan(M,Nt);
S.el    = nan(M,Nt);
S.range = nan(M,Nt);
S.nVis  = zeros(M,Nt);

for m = 1:M
    for k = 1:Nt
        vis = squeeze(G.vis(m,:,k));
        S.nVis(m,k) = sum(vis);             % TODOS los visibles (interferentes)
        vis = vis(:) & satMask;             % solo los ELEGIBLES pueden servir
        if ~any(vis), continue; end
        if strcmpi(cfg.geom.serving,'minRange')
            metric = squeeze(G.range(m,:,k));  metric(~vis) = Inf;  [~,j] = min(metric);
        else
            metric = squeeze(G.el(m,:,k));     metric(~vis) = -Inf; [~,j] = max(metric);
        end
        S.idx(m,k)   = j;
        S.el(m,k)    = G.el(m,j,k);
        S.range(m,k) = G.range(m,j,k);
    end
end
end
