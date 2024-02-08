function X = generateDirectionField(this, x, rootFID, rootVec)
%GENERATEDIRECTIONFIELD This method computes a vector field by parallel
%transport given a trivial connection defined on the dual edges of a mesh
%triangulation. This is based on the implementation in
%'TrivialConnectionsMATLAB' by Keenan Crane 
%
%   INPUT PARAMETERS:
%
%       - x:        #Ex1 trivial connection defined on mesh dual edges
%
%       - rootFID:  The face ID of the root used to populate the output
%                   vector field
%
%       - rootVec:  A 3D unit vector used as the root to populate the
%                   output vector field
%
%   OUTPUT PARAMETRES:
%
%       - X:        #Fx3 output vector field
%
%   by Dillon Cislo 2024/02/07

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if (nargin < 3), rootFID = []; end
if (nargin < 4), rootVec = []; end

validateattributes( x, {'numeric'}, {'vector', 'finite', 'real', ...
    'numel', size(this.E,1)} );
if (size(x,2) ~= 1), x = x.'; end

if isempty(rootFID)
    rootFID = 1; % Choose an arbitrary initial face, if necessary
else
    validateattributes(rootFID, {'numeric'}, {'scalar', 'integer', ...
        'finite', 'real', 'positive', '<=', size(this.F,1)});
end

TR = triangulation(this.F, this.V);
FN = TR.faceNormal;
if isempty(rootVec)
    % Choose an arbitrary initial vector, if necessary
    rootVec = this.V(this.F(rootFID, 1), :) - this.V(this.F(rootFID, 2), :);
    rootVec = rootVec.';
else
    validateattributes(rootVec, {'numeric'}, ...
        {'vector', 'finite', 'real', 'numel', 3});
    if (size(rootVec, 2) ~= 1), rootVec = rootVec.'; end
end
rootVec = rootVec - dot(rootVec, FN(rootFID, :).') .* FN(rootFID, :).';
rootVec = rootVec ./ norm(rootVec);

% Extract the triangle neighbors as a cell array. triNeighbors{f} holds a
% contains the face IDs of all faces that share an edge with face f
triNeighbors = TR.neighbors;
triNeighbors = mat2cell(triNeighbors, ones(size(this.F,1),1), 3);
triNeighbors = cellfun( @(x) x(~isnan(x)), triNeighbors, 'Uni', false);

% Generate a shared edge adjacency matrix. sharedEdgeMatrix(i,j) returns
% the edge ID of the edge shared by face i and face j, and is zero if those
% two faces do not share an edge
sharedEdgeMatrix = TR.edgeAttachments(this.E);
bulkEdges = find(cellfun(@(x) numel(x) > 1, sharedEdgeMatrix, 'Uni', true));
sharedEdgeMatrix = cell2mat(sharedEdgeMatrix(bulkEdges));
sharedEdgeMatrix = sparse(sharedEdgeMatrix, fliplr(sharedEdgeMatrix), ...
    repmat(bulkEdges, 1, 2), size(this.F,1), size(this.F,1));

%--------------------------------------------------------------------------
% Generate Vector Field
%--------------------------------------------------------------------------

X = zeros(3, size(this.F,1));

% Initialize a record of faces that have been visited
visited = false(size(this.F,1), 1);

% Enqueue the root face
X(:, rootFID) = rootVec;
Q = rootFID;
visited(rootFID) = true;

% Transport the initial vector to all other triangles
while ~isempty(Q)
    
    % Dequeue face i
    i = Q(end);
    Q = Q(1:(end-1));
    
    % Visit all neighrbos of face i
    for j = triNeighbors{i}
        
        % Ignore neighbors we've already visited
        if ~visited(j)
            
            % Tranport the tangent vector at face i to face j
            sharedEID = sharedEdgeMatrix(i, j);
            X(:,j) = transportDualTangentVector( ...
                X(:,i), i, j, sharedEID, this.F, this.V.', x);
            
            % Enqueue face j
            Q = [j, Q];
            visited(j) = true;
            
        end
        
    end

end

X = X.';

end

function w = transportDualTangentVector(w0, i, j, eID, F, V, x)
%TRANSPORTDUALTANGENTVECTOR Parallel transport a unit tangent vector from
%face i to face j using the connection specified by x
%
%   INPUT PARAMETERS:
%
%       - w0:       3x1 tangent vector coordinates
%
%       - i:        Index of source triangle
%
%       - j:        Index of destination triangle
%
%       - eID:      Index of shared edge
%
%       - F:        #Fx3 face connectivity list
%
%       - V:        3x#V vertex coordinate list
%
%       - x:        #Ex1 vector giving the deviation from Levi-Civita on
%                   each dual edge
%
%   OUTPUT PARAMETERS:
%
%       - w:        3x1 transported tangent vector coordinates

% Grab vertices a, b, c, and d of the two adjacent
% triangles i and j according to the following labels:
%
%                         b
%                        /|\
%                       / | \
%                      /  |  \
%                     /   |   \
%                    c  i | j  d
%                     \   |   /
%                      \  |  /
%                       \ | /
%                        \|/
%                         a

I = sort(intersect(F(i,:), F(j,:)));
a = V(:, I(1));
b = V(:, I(2));
c = V(:, setdiff(F(i,:), I));
d = V(:, setdiff(F(j,:), I));

% Compute the change of basis between triangles i and j
Ei = orthogonalize( b-a, c-a );
Ej = orthogonalize( b-a, b-d ); % Note the order of entries!

% Compute the in-plane rotation defined by the connection
R = [ cos(x(eID)), -sin(x(eID)), 0; ...
      sin(x(eID)),  cos(x(eID)), 0; ...
               0,             0, 1 ];
           
% Compose these maps to compute the parallel transport
% w = Ej * R * inv(Ei) * w0;
w = Ej * R * (Ei \ w0);

end

function E = orthogonalize(u, v)
%ORTHOGONALIZE Returns a 3x3 orthonormal matrix whose first column is
% parallel to u and whose final column is orthogonal to u and v

   e1 = u;
   e1 = e1/norm(e1);

   e2 = v;
   e2 = e2-(e2'*e1)*e1;
   e2 = e2/norm(e2);

   e3 = cross( e1, e2 );
   e3 = e3/norm(e3);

   E = [ e1, e2, e3 ];
   
end


