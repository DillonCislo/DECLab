%% Discretization Error in the DEC ========================================
%
%   This is a test of how the discretization resolution affects the quality
%   of DEC results. We compare DEC output to analytic results on spherical
%   meshes of increasing resolution
%
%   by Dillon Cislo and Noah P Mitchell 2023
%
%==========================================================================

clear; close all; clc;

% Navigate to tutorial directory
[tutorialDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(tutorialDir)

addpath(genpath('..'));

%% ************************************************************************
% *************************************************************************
%               GENERATE ANAYLTIC RESULTS FOR COMPARISON
% *************************************************************************
% *************************************************************************

syms theta phi x y z
assume( theta, 'real' ); assume( phi, 'real' );
assume( x, 'real'); assume( y, 'real'); assume( z, 'real' );

%==========================================================================
% Generate the Surface Parameters
%==========================================================================

% The equation of the surface
R = [ sin(theta) * cos(phi); sin(theta) * sin(phi); cos(theta) ];

% The tangent vectors
Etheta = [ gradient(R(1), theta); gradient(R(2), theta); gradient(R(3), theta) ];
Ephi = [ gradient(R(1), phi); gradient(R(2), phi); gradient(R(3), phi) ];

% The unit normal vector
% N = cross(Etheta, Ephi);
% N = simplify( N ./ sqrt( sum( N.^2 ) ) );
N = R;

% The metric tensor
g = simplify( [ dot(Etheta, Etheta), dot(Etheta, Ephi); ...
    dot(Ephi, Etheta), dot(Ephi, Ephi) ] );
gInv = simplify(inv(g));

% The dual basis vectors
dtheta = Etheta;
dphi = Ephi ./ sin(theta).^2;

%==========================================================================
% Generate a Scalar Field on the Surface
%==========================================================================

%--------------------------------------------------------------------------
% Enter your favorite scalar field in Cartesian or spherical coordinates

S = cos(4*phi).*cos(theta).*sin(theta).^4-cos(theta).^2;

%--------------------------------------------------------------------------

% Transform to spherical coordinates if necessary
S = simplify(subs( S, [x y z], [R(1) R(2) R(3)] ));

% Calculate the gradient of the scalar field
gradS = simplify( gradient(S, theta) * Etheta + ...
    ( gradient(S, phi) / sin(theta) ) * ( Ephi / sin(theta) ) );

% Calculate the Laplacian of the scalar field
lapS = simplify( ...
    gradient( sin(theta) * gradient(S, theta), theta ) / sin(theta) + ...
    gradient( gradient( S, phi ), phi ) / sin(theta)^2 );

%==========================================================================
% Generate a Vector Field on the Surface
%==========================================================================

%--------------------------------------------------------------------------
% Enter your favorite vector field in Cartesian or spherical coordinates

U = [ x * z * ( z^2 - 1/4 ) - y; ...
    y * z * ( z^2 - 1/4 ) + x; ...
    -( x^2 + y^2 ) * ( z^2 - 1/4 ) ];

%--------------------------------------------------------------------------

% Transform to spherical coordinates if necessary
U = simplify(subs( U, [x, y, z], [R(1) R(2) R(3)] ));

% Project onto the tangent space of the surface if necessary
U = U - dot(U, N) * N;

% Calculate the components of U in the (theta, phi)-basis
Utheta = simplify( dot(U, Etheta) / dot(Etheta, Etheta) );
Uphi = simplify( dot(U, Ephi) ./ dot(Ephi, Ephi) );

% Calculate the components of U in the dual basis
uTheta = Utheta;
uPhi = sin(theta)^2 * Uphi;

% Calculate the divergence of the vector field
divU = simplify( gradient( sin(theta) * uTheta, theta ) / sin(theta) + ...
    gradient( uPhi, phi ) / sin(theta)^2 );

% Calculate the 'curl' of the vector field
curlU = simplify(...
    (gradient( uPhi, theta ) - gradient( uTheta, phi )) / sin(theta) );

% Calculate the Laplacian of the vector field -----------------------------

 % Exterior derivative of the associated 1-form field
dU = (gradient(uPhi, theta)-gradient(uTheta, phi)) / sqrt(det(g));

% d * d * U
lapU1 = ...
    gradient(sqrt(det(g)) * (uTheta * gInv(1,1) + uPhi * gInv(2,1)), theta) + ...
    gradient(sqrt(det(g)) * (uTheta * gInv(1,2) + uPhi * gInv(2,2)), phi);
lapU1 = simplify(lapU1 / sqrt(det(g)));
lapU1 = gradient(lapU1, theta) * dtheta + gradient(lapU1, phi) * dphi;
lapU1 = simplify(lapU1);

%  * d * d U
lapU2 = ...
    (gradient(dU, theta) * gInv(1,2) + gradient(dU, phi) * gInv(2,2)) * dtheta - ...
    (gradient(dU, theta) * gInv(1,1) + gradient(dU, phi) * gInv(2,1)) * dphi;
lapU2 = simplify(-sqrt(det(g)) * lapU2);

lapU = simplify(lapU1 + lapU2);

% Covariant 1-form components in dual basis
lapUTheta = simplify( dot(lapU, Etheta) ); 
lapUPhi = simplify( dot(lapU, Ephi) );

% Convert 1-form field to tangent vector field
lapU = (gInv(1,1) * lapUTheta + gInv(1,2) * lapUPhi) * Etheta + ...
    (gInv(2,1) * lapUTheta + gInv(2,2) * lapUPhi) * Ephi;


%% Clear Extraneous Variables =============================================

clear R Etheta Ephi N g dtheta dphi Utheta Uphi uTheta uPhi
clear dU lapU1 lapU2 lapUTheta lapUPhi

%% ************************************************************************
% *************************************************************************
%                       TEST DIFFERENTIAL OPERATORS
% *************************************************************************
% *************************************************************************

%% Calculate the Gradient of a Scalar Field ===============================
% The gradient calculated by DEC exactly matches the classical Finite
% Element Method (FEM) gradient (see 'grad.m' from 'gptoolbox')

close all; clc;

triVertexNum = 2.^(5:14);
numFaces = zeros(size(triVertexNum));
allRelErr = cell(numel(triVertexNum), 1);

for i = 1:numel(triVertexNum)
    
    fprintf('Processing surface %d/%d... ', i, numel(triVertexNum));
    
    % Produce the new triangulation
    sphereTri = sphereTriangulationVogel(triVertexNum(i));
    F = sphereTri.ConnectivityList;
    V = sphereTri.Points;
    
    numFaces(i) = size(F,1); % Number of faces in the current triangulation
    
    % Evaluate analytic values on current triangulation -------------------
    
    % (theta, phi) for each vertex
    NTheta = acos(V(:,3));
    NPhi = atan2(V(:,2), V(:,1));
    
    anaS = matlabFunction(S, 'Vars', {theta, phi});
    anaS = anaS(NTheta, NPhi);
    
    anaGradS = matlabFunction(gradS.', 'Vars', {theta, phi});
    anaGradS = anaGradS(NTheta, NPhi);
    
    % Average vector fields onto faces
    anaGradS = cat(3, anaGradS(F(:,1), :), anaGradS(F(:,2), :), ...
        anaGradS(F(:,3), :) );
    anaGradS = mean( anaGradS, 3 );
    
    % Perform DEC analysis ------------------------------------------------
    
    DEC = DiscreteExteriorCalculus( F, V );
    NGradS = DEC.gradient(anaS);
    
    % The relative error
    relErr = anaGradS - NGradS;
    relErr = sqrt(sum(relErr.^2, 2)) ./ sqrt(sum(anaGradS.^2, 2));
    
    % Account for division by 0
    relErr(isinf(relErr)) = 0;
    relErr(isnan(relErr)) = 0;
    
    allRelErr{i} = relErr;
    
    fprintf('Done\n');
    
end

allMedErr = cellfun(@median, allRelErr, 'Uni', true);
plot(numFaces, allMedErr, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'b');
set(gca, 'XScale', 'log');
xlabel('Number of Mesh Faces');
ylabel('Median Fractional Error in \nabla S');

clear i sphereTri F V NTheta NPhi anaS anaGradS DEC NGradS relErr
clear triVertexNum numFaces allRelErr allMedErr allRMSErr

%% Calculate the Laplacian of a Scalar Field ==============================
% Without area weight re-normalization, the DEC Laplacian is identical to
% the cotangent-weight classical FEM Laplacian (see 'cotmatrix.m' from
% 'gptoolbox').  With area-weight renormalization the DEC Laplcian matches
% Eqn. (3.11) from "Polygon Mesh Processing" by Botsch et al.  The latter
% formalism seems to be superior in so far as matching analytic results

close all; clc;

triVertexNum = 2.^(5:14);
numFaces = zeros(size(triVertexNum));
allRelErr = cell(numel(triVertexNum), 1);
normalizeAreas = true;

for i = 1:numel(triVertexNum)
    
    fprintf('Processing surface %d/%d... ', i, numel(triVertexNum));
    
    % Produce the new triangulation
    sphereTri = sphereTriangulationVogel(triVertexNum(i));
    F = sphereTri.ConnectivityList;
    V = sphereTri.Points;
    
    numFaces(i) = size(F,1); % Number of faces in the current triangulation
    
    % Evaluate analytic values on current triangulation -------------------
    
    % (theta, phi) for each vertex
    NTheta = acos(V(:,3));
    NPhi = atan2(V(:,2), V(:,1));
    
    anaS = matlabFunction(S, 'Vars', {theta, phi});
    anaS = anaS(NTheta, NPhi);
    
    anaLapS = matlabFunction(lapS, 'Vars', {theta, phi});
    anaLapS = anaLapS(NTheta, NPhi);
    
    % Perform DEC analysis ------------------------------------------------
    
    DEC = DiscreteExteriorCalculus( F, V );
    NLapS = DEC.laplacian(anaS, normalizeAreas);
    
    % The relative error
    relErr = anaLapS - NLapS;
    relErr = sqrt(sum(relErr.^2, 2)) ./ sqrt(sum(anaLapS.^2, 2));
    
    % Account for division by 0
    relErr(isinf(relErr)) = 0;
    relErr(isnan(relErr)) = 0;
    
    allRelErr{i} = relErr;
    
    fprintf('Done\n');
    
end

allMedErr = cellfun(@median, allRelErr, 'Uni', true);
plot(numFaces, allMedErr, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'b');
set(gca, 'XScale', 'log');
xlabel('Number of Mesh Faces');
ylabel('Median Fractional Error in \nabla^2 S');

clear i sphereTri F V NTheta NPhi anaS anaLapS DEC NLapS relErr
clear triVertexNum numFaces allRelErr allMedErr allRMSErr  normalizeAreas

%% Calculate the Divergence of a Tangent Vector Field =====================
% The divergence calculated by the DEC does NOT match the classical FEM
% divergence (or at least the FEM divergence implemented in 'div.m'
% in 'gptoolbox' - it might match an FEM divergence with different
% weights), HOWEVER it does seem to perform significantly better than the
% FEM divergence in so far as matching analytic results!

close all; clc;

triVertexNum = 2.^(5:14);
numFaces = zeros(size(triVertexNum));
allRelErr = cell(numel(triVertexNum), 1);

for i = 1:numel(triVertexNum)
    
    fprintf('Processing surface %d/%d... ', i, numel(triVertexNum));
    
    % Produce the new triangulation
    sphereTri = sphereTriangulationVogel(triVertexNum(i));
    F = sphereTri.ConnectivityList;
    V = sphereTri.Points;
    
    numFaces(i) = size(F,1); % Number of faces in the current triangulation
    
    % Evaluate analytic values on current triangulation -------------------
    
    % (theta, phi) for each vertex
    NTheta = acos(V(:,3));
    NPhi = atan2(V(:,2), V(:,1));
    
    anaU = matlabFunction(U.', 'Vars', {theta, phi});
    anaU = anaU(NTheta, NPhi);
    
    anaDivU = matlabFunction(divU, 'Vars', {theta, phi});
    anaDivU = anaDivU(NTheta, NPhi);
    
    % Perform DEC analysis ------------------------------------------------
    
    DEC = DiscreteExteriorCalculus( F, V );
    NDivU = DEC.divergence(anaU);
    
    % The relative error
    relErr = anaDivU - NDivU;
    relErr = sqrt(sum(relErr.^2, 2)) ./ sqrt(sum(anaDivU.^2, 2));
    
    % Account for division by 0
    relErr(isinf(relErr)) = 0;
    relErr(isnan(relErr)) = 0;
    
    allRelErr{i} = relErr;
    
    fprintf('Done\n');
    
end

allMedErr = cellfun(@median, allRelErr, 'Uni', true);
plot(numFaces, allMedErr, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'b');
set(gca, 'XScale', 'log');
xlabel('Number of Mesh Faces');
ylabel('Median Fractional Error in \nabla \cdot U');

clear i sphereTri F V NTheta NPhi anaU anaDivU DEC NDivU relErr
clear triVertexNum numFaces allRelErr allMedErr allRMSErr

%% Calculate the "Curl" of a Tangent Vector Field =========================
% There is no correspondence to the "curl" of a tangent vector field in the
% classical FEM formalism that manages to treat the "curl" as a primal
% 2-form on mesh facets.  You might be able to compare the results to some
% FEM-style measurement of circulation on a hinge or vertex-stencil, but
% these results seem convincing enough me without such a comparison.

close all; clc;

triVertexNum = 2.^(5:14);
numFaces = zeros(size(triVertexNum));
allRelErr = cell(numel(triVertexNum), 1);

for i = 1:numel(triVertexNum)
    
    fprintf('Processing surface %d/%d... ', i, numel(triVertexNum));
    
    % Produce the new triangulation
    sphereTri = sphereTriangulationVogel(triVertexNum(i));
    F = sphereTri.ConnectivityList;
    V = sphereTri.Points;
    
    numFaces(i) = size(F,1); % Number of faces in the current triangulation
    
    % Evaluate analytic values on current triangulation -------------------
    
    % (theta, phi) for each vertex
    NTheta = acos(V(:,3));
    NPhi = atan2(V(:,2), V(:,1));
    
    anaU = matlabFunction(U.', 'Vars', {theta, phi});
    anaU = anaU(NTheta, NPhi);
    
    % Average vector fields onto faces
    anaU = cat(3, anaU(F(:,1), :), anaU(F(:,2), :), ...
        anaU(F(:,3), :) );
    anaU = mean( anaU, 3 );
    
    anaCurlU = matlabFunction(curlU, 'Vars', {theta, phi});
    anaCurlU = anaCurlU(NTheta, NPhi);
    
    % Perform DEC analysis ------------------------------------------------
    
    DEC = DiscreteExteriorCalculus( F, V );
    NCurlU = DEC.curl(anaU, 'dual');
    
    % The relative error
    relErr = anaCurlU - NCurlU;
    relErr = sqrt(sum(relErr.^2, 2)) ./ sqrt(sum(anaCurlU.^2, 2));
    
    % Account for division by 0
    relErr(isinf(relErr)) = 0;
    relErr(isnan(relErr)) = 0;
    
    allRelErr{i} = relErr;
    
    fprintf('Done\n');
    
end

allMedErr = cellfun(@median, allRelErr, 'Uni', true);
plot(numFaces, allMedErr, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'b');
set(gca, 'XScale', 'log');
xlabel('Number of Mesh Faces');
ylabel('Median Fractional Error in \nabla \times U');

clear i sphereTri F V NTheta NPhi anaU anaCurlU DEC NCurlU relErr
clear triVertexNum numFaces allRelErr allMedErr allRMSErr

%% Calculate the Laplacian of a Tangent Vector Field ======================
% NOTE: This functionality is somewhat experimental. 
%
% We define the Laplace-de Rham operator for k-forms as (*d*d + d*d*) u.
% For some reason the pathway (d * d * u) works for primal 1-forms, but
% (* d * d u) does not. Conversely, (* d * d u) works for dual 1-forms, but
% (d * d * u) does not. It is not immediately clear why, but it seems to be
% linked to the application of the 'dd0' operator that transforms dual
% 0-forms to dual 1-forms by exterior differentiation. Again, this is
% strange since this operator is confirmed to have the same form as
% multiple other sources/DEC packages.
%
% The Laplacian for a dual vector field U is calculated by
% converting it to a primal 1-form up and calculating (d * d * up), then
% converting U to a dual 1-form ud and calculating (* d * d ud), then
% converting both quantities back to dual vector fields and adding them
% together
%
% The results of this calculation are the least accurate relative to
% analytic results of any of our DEC operators (~2% median error at ~30k
% faces) and appears to be highly mesh dependent. More work needs to be
% done to diagnose why. Use at your own risk.

close all; clc;

triVertexNum = 2.^(5:14);
numFaces = zeros(size(triVertexNum));
allRelErr = cell(numel(triVertexNum), 1);

for i = 1:numel(triVertexNum)
    
    fprintf('Processing surface %d/%d... ', i, numel(triVertexNum));
    
    % Produce the new triangulation
    sphereTri = sphereTriangulationVogel(triVertexNum(i));
    F = sphereTri.ConnectivityList;
    V = sphereTri.Points;
    FN = sphereTri.faceNormal;
    
    numFaces(i) = size(F,1); % Number of faces in the current triangulation
    
    % Evaluate analytic values on current triangulation -------------------
    
    % (theta, phi) for each vertex
    NTheta = acos(V(:,3));
    NPhi = atan2(V(:,2), V(:,1));
    
    anaU = matlabFunction(U.', 'Vars', {theta, phi});
    anaU = anaU(NTheta, NPhi);
    
    % Average vector fields onto faces
    anaU = cat(3, anaU(F(:,1), :), anaU(F(:,2), :), ...
        anaU(F(:,3), :) );
    anaU = mean( anaU, 3 );
    anaU = anaU - repmat(dot(anaU, FN, 2), 1, 3) .* FN;
    
    anaLapU = matlabFunction(lapU.', 'Vars', {theta, phi});
    anaLapU = anaLapU(NTheta, NPhi);
    
    % Average vector fields onto faces
    anaLapU = cat(3, anaLapU(F(:,1), :), anaLapU(F(:,2), :), ...
        anaLapU(F(:,3), :) );
    anaLapU = mean( anaLapU, 3 );
    anaLapU = anaLapU - repmat(dot(anaLapU, FN, 2), 1, 3) .* FN;
    
    % Perform DEC analysis ------------------------------------------------
    
    DEC = DiscreteExteriorCalculus( F, V );
    NLapU = DEC.laplacian(anaU);
    
    % The relative error
    relErr = anaLapU - NLapU;
    relErr = sqrt(sum(relErr.^2, 2)) ./ sqrt(sum(anaLapU.^2, 2));
    
    % Account for division by 0
    relErr(isinf(relErr)) = 0;
    relErr(isnan(relErr)) = 0;
    
    allRelErr{i} = relErr;
    
    fprintf('Done\n');
    
end

allMedErr = cellfun(@median, allRelErr, 'Uni', true);
plot(numFaces, allMedErr, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'b');
set(gca, 'XScale', 'log');
xlabel('Number of Mesh Faces');
ylabel('Median Fractional Error in \nabla^2 U');

clear i sphereTri F V NTheta NPhi anaU anaLapU DEC NLapU relErr
clear triVertexNum numFaces allRelErr allMedErr allRMSErr
