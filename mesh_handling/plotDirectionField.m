function fig = plotDirectionField(F, V, X, defectIDx, ...
    plotOptions, scatterOptions)
%PLOTDIRECTIONFIELD Plots a tangent direction field on a mesh. Also plots
%singular vertices. Primarily intended to quickly view the results of the
%'computeTrivialConnection' method
%
%   INPUT PARAMETERS:
%
%       - F:            #Fx3 face connectivity list
%
%       - V:            #Vx3 vertex coordinate list
%
%       - X:            #Fx3 dual tangent vector field specifying a unit
%                       vector in each face
%
%       - defectIDx:    #Dx1 vector of vertex indices denoting the
%                       locations of the user specified singularites
%
%       - plotOptions:  A cell array holding the (name, value) pair options
%                       used to generate the direction field plot
%
%       - scatterOptions:   A cell array holding the (name, value) pair
%                           options used to plot the singularities
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

validateattributes(V, {'numeric'}, {'2d', 'ncols', 3, 'finite', 'real'});
validateattributes(F, {'numeric'}, {'2d', 'ncols', 3, 'finite', ...
    'integer', 'positive', 'real', '<=', size(V,1)});

validateattributes(X, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real', 'nrows', size(F,1)});
X = X ./ sqrt(sum(X.^2, 2)); % Normalize just in case

if ~isempty(defectIDx)
    validateattributes(defectIDx, {'numeric'}, {'vector', 'integer', ...
        'positive', 'finite', 'real', '<=', size(V,1)} );
    if (size(defectIDx,2) ~= 1), defectIDx = defectIDx.'; end
end

if (nargin < 5)
    plotOptions = {'Color', [0 0 1], 'LineWidth', 2};
end

if (nargin < 6)
    scatterOptions = {'Color', 'r', 'MarkerSize', 40};
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
a = COM + R .* X + smallDist * FN;
b = COM - R .* X + smallDist * FN;
% c = nan(size(a));

% plotPoints = [a; b];
% plotPoints = nan(3*size(F,1), 3);
% plotPoints(1:3:end) = COM + R .* X + smallDist * FN;
% plotPoints(2:3:end) = COM - R .* X + smallDist * FN;

%--------------------------------------------------------------------------
% Generate Visualization
%--------------------------------------------------------------------------

fig = figure('Color', 'w');

trisurf(triangulation(F, V), 'FaceColor', 0.9 * ones(1,3));
hold on
plot3([a(:,1).'; b(:,1).'], [a(:,2).'; b(:,2).'], [a(:,3).'; b(:,3).'], ...
    'Color', [0 0 1], 'LineWidth', 2);
scatter3(V(defectIDx, 1), V(defectIDx, 2), V(defectIDx, 3), 40, 'filled', 'r');
hold off

axis equal off
cameratoolbar('SetMode', 'orbit');



end

