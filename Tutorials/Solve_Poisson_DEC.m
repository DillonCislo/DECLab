%% Solving the Poisson Equation with the DEC ==============================
% This is a script to show how to solve the Poisson equation with Dirichlet
% boundary conditions on an unstructured 2D domain using the Discrete
% Exterior Calculus (DEC) formalism
%
% by Dillon Cislo 2020/11/15
%==========================================================================
clear; close all; clc;

%--------------------------------------------------------------------------
% Construct a Basic 2D Triangulation
%--------------------------------------------------------------------------

% A triangulation of the unit disk ----------------------------------------
TR = diskTriangulation(30);

F = TR.ConnectivityList; % Face connectivity list
V = TR.Points; % Vertex coordinate list

% A rectilinear grid ------------------------------------------------------

% % Create basic grid points in 2D
% L = 2; % The length of the domain
% W = 2; % The width of the domain
% NL = 50; % The number of vertices along the length of the cylinder
% NW = round( W * NL ./ L); % The number of vertices along the width
% [X, Y] = meshgrid( linspace(0, L, NL), linspace(0, W, NW) );
%
% % Center on (0,0)
% X = X-L/2;
% Y = Y-W/2;
%
% % Triangulate grid points in 2D
% TR = delaunayTriangulation([X(:), Y(:)]);
%
% F = bfs_orient(TR.ConnectivityList);
% V = TR.Points;
%
% clear L NL NW Y X

%--------------------------------------------------------------------------
% Re-Format the Vertex Coordinate List so that Boundary Vertices are at the
% End of the List
%--------------------------------------------------------------------------

bdyIDx = unique(freeBoundary(TR));

newVIDx = (1:size(V,1)).';
newVIDx(bdyIDx) = [];
newVIDx = [newVIDx; bdyIDx];

V = V(newVIDx, :);
F = changem( F, (1:size(V,1)).', newVIDx );

TR = triangulation(F, V);
E = TR.edges; % Edge Connectivity List

bdyIDx = unique(freeBoundary(TR));

clear newVIDx

%% Generate Analytic Results ==============================================

syms x y
assume( x, 'real' );
assume( y, 'real' );

% Construct the function that will be solved for
u = y * (1 - y) * x^3;

% Calculate the Laplacian of the function
g = gradient(gradient(u,x),x) + gradient(gradient(u,y),y);

%% Convert Symbolic Quantities to Numerical Quantities ====================

fprintf('Substituting numerical values for symbolic variables... ');

X = V(:,1); Y = V(:,2);

U = double(vpa(subs(u, {x,y}, {X,Y})));
G = double(vpa(subs(g, {x,y}, {X,Y})));

fprintf('Done\n');

%% Solve the Poisson Equation =============================================

% Construct Differential Operators ----------------------------------------

% A DEC object for the current mesh
DEC = DiscreteExteriorCalculus( F, [ V, zeros(size(V,1), 1) ] );

% The full mesh vertex 'mass' operator
M = DEC.hd0;

% The full mesh (unweighted) Laplacian matrix
% ( The weighted Laplacian L = M^(-1) * C )
C = DEC.dd1 * DEC.hd1 * DEC.d0;

% NOTE: IT IS REQUIRED THAT ALL THE BOUNDARY VERTICES BE AT THE END OF
% THE VERTEX COORDINATE LIST (THIS SHOULD BE DONE DURING MESH
% CONSTRUCTION)

% The 'mass' operator for interior vertices
MII = M(1:(min(bdyIDx)-1), 1:(min(bdyIDx)-1));

% The (unweighted) Laplacian for interior vertices
CII = C(1:(min(bdyIDx)-1), 1:(min(bdyIDx)-1));

% The mixed (unweighted) Laplacian
CIB = C(1:(min(bdyIDx)-1), min(bdyIDx):end);

% Construct the Poisson Kernel --------------------------------------------

% The kernel for interior vertices
GII = G((1:(min(bdyIDx)-1)).');

% Set Dirichlet Boundary Conditions ---------------------------------------

% The function values for boundary vertices
UB = U(bdyIDx);

%--------------------------------------------------------------------------
% Solve the Poisson Problem
%--------------------------------------------------------------------------

% Calculate the velocity components
calcU = CII \ (MII * GII - CIB * UB);

% Calculate solution residuals (checks if any solution was found - not that
% the solution is correct)
solveRes = ( CII * calcU - (MII * GII - CIB * UB) );

% Include known boundary velocities
calcU = [ calcU; UB ];

fprintf('Maximum Solution Residual = %0.5e\n', max(abs(solveRes)));

%% Check Results ==========================================================
close all; clc;

poissonErr = abs(U - calcU);

fprintf('Maximum Error = %0.5e\n', max(poissonErr));
fprintf('RMS Error = %0.5e\n', sqrt(mean(poissonErr.^2)));
fprintf('Median Error = %0.5e\n', median(poissonErr));

figure

subplot(1,3,1)

patch('Faces', F, 'Vertices', V, 'FaceVertexCData', U, ...
    'FaceColor', 'interp', 'EdgeColor', 'none');

axis equal
colorbar

title('The True Function');

subplot(1,3,2)

patch('Faces', F, 'Vertices', V, 'FaceVertexCData', calcU, ...
    'FaceColor', 'interp', 'EdgeColor', 'none');

axis equal
colorbar
title('The Calculated Solution');

subplot(1,3,3)

patch('Faces', F, 'Vertices', V, 'FaceVertexCData', poissonErr, ...
    'FaceColor', 'interp', 'EdgeColor', 'none');

axis equal
colorbar
title('The Error');


