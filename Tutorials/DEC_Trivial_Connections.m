%% DEC Trivial Connections Test ===========================================
%
%   This script shows how to use 'Discrete Exterior Calculus' class to
%   compute trivial connections and direction fields on genus 0 surfaces
%
%   by Dillon Cislo 2024
%
%==========================================================================

[tutorialDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(tutorialDir)
addpath('..');
addpath('../mesh_handling');
addpath('../PlottingFunctions');


%% Spherical Mesh Example =================================================
clear; close all; clc;

[tutorialDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(tutorialDir);

% Generate a mesh of the unit sphere
sphereTri = sphereTriangulationVogel(750);
F = sphereTri.ConnectivityList; V = sphereTri.Points;
DEC = DiscreteExteriorCalculus(F, V);

% Specify desired singularities. NOTE: Defect charge *must* sum to 2
% defectIDx = knnsearch(V, [0 0 1]);
% defectCharge = 2;

% defectIDx = knnsearch(V, [0 0 1; 0 0 -1]);
% defectCharge = [1, 1];

% defectIDx = knnsearch(V, [0 0 1; sqrt(2)/3 0 -0.5; -sqrt(2)/3 0 -0.5]);
% defectCharge = [1, -1, 2];

defectIDx = randsample(1:size(V,1), 3);
% defectIDx = knnsearch(V, [0 0 1; sqrt(2)/3 0 -0.5; -sqrt(2)/3 0 -0.5]);
defectCharge = [0.5, -0.5, 2];

% Compute the trivial connection
x = DEC.computeTrivialConnection(defectIDx, defectCharge);

% Generate a corresponding vector field
U = DEC.generateDirectionField(x);

% Plot results
plotDirectionField(F, V, U, defectIDx);

%% Disk Mesh Example ======================================================
clear; close all; clc;

[tutorialDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(tutorialDir);

diskTri = diskTriangulation(20);
F = diskTri.ConnectivityList; V = diskTri.Points;
V = [V, zeros(size(V,1), 1)];
% V = [V, exp(-sum(V.^2, 2))];
DEC = DiscreteExteriorCalculus(F, V);

% Specify desired singularities. There are no requirements on the total sum
% of the defect charge for a disk
defectIDx = knnsearch(V, [0 0 1]);
defectCharge = 1;

% defectIDx = knnsearch(V, [0.5 0 0; -0.5 0 0]);
% defectCharge = [+1/2, -1/2];

% Compute the trivial connection
x = DEC.computeTrivialConnection(defectIDx, defectCharge);

% Generate a corresponding vector field
U = DEC.generateDirectionField(x);
% U = DEC.helmholtzHodgeDecomposition(U, 1e-8);

% Plot results
plotDirectionField(F, V(:, 1:2), U(:, 1:2), defectIDx);
% view([0, 90]);

