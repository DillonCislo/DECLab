function fig = plotDirectionField(F, V, X, defectIDx, ...
    patchOptions, plotOptions, scatterOptions, visFIDx)
%PLOTDIRECTIONFIELD Plots a tangent direction field on a mesh. Also plots
%singular vertices. Primarily intended to quickly view the results of the
%'computeTrivialConnection' method
%
%   INPUT PARAMETERS:
%
%       - F:            #Fx3 face connectivity list
%
%       - V:            #VxD vertex coordinate list
%
%       - X:            #Fx3 dual tangent vector field specifying a unit
%                       vector in each face
%
%       - defectIDx:    #Dx1 vector of vertex indices denoting the
%                       locations of the user specified singularites
%
%       - patchOptions:     A cell array holding the (name, value)-pair
%                           options used to generate the surface plot
%
%       - plotOptions:      A cell array holding the (name, value)-pair
%                           options used to generate the direction field
%                           plot
%
%       - scatterOptions:   A cell array holding the (name, value)-pair
%                           options used to plot the singularities
%
%       - visFIDx:          #VFx1 vector of face IDs on which to view the
%                           direction field (default is all faces)
%
%   OUTPUT PARAMETERS:
%
%       - fig:          The handle to the generated figure
%
% by Dillon Cislo 2024/02/07

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), error('Please supply vertex coordinate list'); end
if (nargin < 3), error('Please supply direction field'); end
if (nargin < 4), defectIDx = []; end

is2D = false;
validateattributes(V, {'numeric'}, {'2d', 'finite', 'real'});
if (size(V,2) == 2)
    is2D = true;
    V = [V, zeros(size(V,1),1)];
elseif (size(V,2) ~= 3)
    error('Vertex coordinates must be 2D or 3D');
end

validateattributes(F, {'numeric'}, {'2d', 'ncols', 3, 'finite', ...
    'integer', 'positive', 'real', '<=', size(V,1)});

validateattributes(X, {'numeric'}, ...
    {'2d', 'finite', 'real', 'nrows', size(F,1)});
X = X ./ sqrt(sum(X.^2, 2)); % Normalize just in case
if (size(X,2) == 2)
    assert(is2D, 'Dimensions of direction field must match vertices');
    X = [X, zeros(size(X,1),1)];
elseif (size(X,2) ~= 3)
    error('Direction field must be 2D or 3D');
end
   

if ~isempty(defectIDx)
    validateattributes(defectIDx, {'numeric'}, {'vector', 'integer', ...
        'positive', 'finite', 'real', '<=', size(V,1)} );
    if (size(defectIDx,2) ~= 1), defectIDx = defectIDx.'; end
end

if (nargin < 5)
    patchOptions = {'FaceColor',  0.9 * ones(1,3)};
elseif isempty(patchOptions)
    patchOptions = {'FaceColor',  0.9 * ones(1,3)};
else
    assert(iscell(patchOptions), 'Patch options must be a cell array');
end

if (nargin < 6)
    plotOptions = {'Color', [0 0 1], 'LineWidth', 2};
elseif isempty(plotOptions)
    plotOptions = {'Color', [0 0 1], 'LineWidth', 2};
else
    assert(iscell(plotOptions), 'Line plot options must be a cell array');
end

if (nargin < 7)
    scatterOptions = {'MarkerFaceColor', 'r', 'SizeData', 40};
elseif isempty(scatterOptions)
    scatterOptions = {'MarkerFaceColor', 'r', 'SizeData', 40};
else
    assert(iscell(scatterOptions), ...
        'Scatter plot options must be a cell array');
end

if (nargin < 8)
    visFIDx = (1:size(F,1)).';
elseif isempty(visFIDx)
    visFIDx = (1:size(F,1)).';
else
    validateattributes(visFIDx, {'numeric'}, {'vector', ...
        'integer', 'positive', 'finite', 'real', '<=', size(F,1)});
    if (size(visFIDx,2) ~= 1), visFIDx = visFIDx.'; end
end

TR = triangulation(F, V);
FN = TR.faceNormal;

% Compute face barycenters
COM = cat(3, V(F(:,1), :), V(F(:,2), :), V(F(:,3), :));
COM = mean(COM, 3);

% Compute the lengths of the edges in each face
l1 = sqrt(sum((V(F(:,3), :) - V(F(:,2), :)).^2, 2));
l2 = sqrt(sum((V(F(:,1), :) - V(F(:,3), :)).^2, 2));
l3 = sqrt(sum((V(F(:,2), :) - V(F(:,1), :)).^2, 2));

% Compute the inradius of each face
R = sqrt( ((l1+l2-l3) .* (l3+l1-l2) .* (l2+l3-l1)) ./ (l1+l2+l3) ) ./ 2;

% Reformat this data to feed into the plotting functions
smallDist = 1e-3;
if is2D
    a = COM + R .* X;
    b = COM - R .* X;
else
    a = COM + R .* X + smallDist * FN;
    b = COM - R .* X + smallDist * FN;
end

%--------------------------------------------------------------------------
% Generate Visualization
%--------------------------------------------------------------------------

if isempty(findall(0, 'Type', 'figure'))
    fig = figure('Color', 'w');
else
    fig = gcf;
end

if is2D
    
    patch('Faces', F, 'Vertices', V, patchOptions{:});
    
    hold on
    
    plot([a(visFIDx, 1).'; b(visFIDx, 1).'], ...
        [a(visFIDx, 2).'; b(visFIDx, 2).'], plotOptions{:} );
    
    scatter(V(defectIDx, 1), V(defectIDx, 2), scatterOptions{:} );
    
    hold off
    
    view(0, 90);
    
else
    
    trisurf(triangulation(F,V), patchOptions{:});
    
    hold on
    
    plot3([a(visFIDx, 1).'; b(visFIDx, 1).'], ...
        [a(visFIDx, 2).'; b(visFIDx, 2).'], ...
        [a(visFIDx, 3).'; b(visFIDx, 3).'],  plotOptions{:} );
    
    scatter3(V(defectIDx, 1), V(defectIDx, 2), V(defectIDx, 3), ...
        scatterOptions{:} );
    
    hold off
    
    cameratoolbar('SetMode', 'orbit');
    
end

axis equal off

end

