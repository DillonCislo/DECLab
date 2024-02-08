function x = computeTrivialConnection(this, defectIDx, defectCharge)
%COMPUTETRIVIALCONNECTION Computes a set of connection angles x associated
%with each dual edge of a mesh that describe a trivial connection with
%the specified singularities. This angles are expressed relative to the
%discrte Levi-Civita connection, i.e. ``translation'' across the shared
%edge between two faces. This is an implementation of algorithm in
%``Trivial connections on discrete surfaces'' by Crane et al. (2010),
%simplified for simply connected disks and spheres (see the note posted at
%https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/ for more
%details)
%
%   INPUT PARAMETERS:
%
%       - defectIDx:    #Dx1 vector of vertex indices denoting the
%                       locations of the user specified singularites
%
%       - defectCharge: #Dx1 vector of the defect charge. On a sphere,
%                       these must sum to 2
%
%   OUTPUT PARAMETERS:
%
%       - x:    #Ex1 vector of angles encoding deviations from the
%               Levi-Civita connection
%
%   by Dillon Cislo 2024/02/07

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
validateattributes(defectIDx, {'numeric'}, {'vector', 'integer', ...
    'positive', 'finite', 'real', '<=', size(this.V,1)} );
if (size(defectIDx,2) ~= 1), defectIDx = defectIDx.'; end

validateattributes(defectCharge, {'numeric'}, ...
    {'vector', 'finite', 'real', 'numel', numel(defectIDx)});
if (size(defectCharge,2) ~= 1), defectCharge = defectCharge.'; end

eulerChi = size(this.F,1) + size(this.V,1) - size(this.E,1);
if (numel(this.allBdys) == 0)
    
    assert( eulerChi == 2, ['This method only works for topological ' ...
        'disks and spheres'] );
    assert( sum(defectCharge) == 2, ['Defect charge must sum to 2 ' ...
        'for topological spheres' ] );
    
elseif (numel(this.allBdys) == 1)
    
    assert( eulerChi == 1, ['This method only works for topological ' ...
        'disks and spheres'] );
    
else
    
    error('This method only works for topological disks and spheres');
    
end

%--------------------------------------------------------------------------
% Compute Connection
%--------------------------------------------------------------------------

% Construct the POSITIVE definite Laplace-Beltrami operator (note the sign
% here!)
L = this.d0.' * this.hd1 * this.d0;

% The right hand side of the system is the Riemannian holonomy around
% vertices (i.e., Gaussian curvature) minus the target curvature at
% singular vertices
b = this.discreteGaussianCurvature();
b(defectIDx) = b(defectIDx) - 2 * pi * defectCharge;

% Compute this minimum norm solution to the problem Lu = -b
u = lsqminnorm(L, -b);

% Compute the trivial connection
x = this.hd1 * this.d0 * u;

end