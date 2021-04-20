function constructDECOperators(this)
%CONSTRUCTDECOPERATORS Construct the various discrete exterior derivative
%operators and discrete Hodge dual operators necessary to utilize the
%Discrete Exterior Calculus on a flat/curved mesh triangulation.
%
%   INTENDED FOR INTERNAL USE WITH 'DiscreteExteriorCalculus.m'
%   WARNING: NO INPUT CHECKS ARE PERFORMED
%
%   by Dillon Cislo 01/29/2020

%==========================================================================
% Process Input Triangulation
%==========================================================================
F = this.F; % Face connectivity list
V = this.V; % Vertex coordinate list
E = this.E; % Edge connectivity list

% A MATLAB-style representation of the triangulation
TR = triangulation(F, V);

numE = size(E,1); % The number of edges
numF = size(F,1); % The number of faces
numV = size(V,1); % The number of vertices

%--------------------------------------------------------------------------
% Construct Topological Structure Tools
%--------------------------------------------------------------------------

% #Ex2 array of fIDs of the faces attached to a particular edge.
% If an edge is a border edge (i.e., only attached to a single face), then
% that fID is listed twice for dimensional consistency
resizeCell = @(x) repmat( x, 1, 1+mod(numel(x),2) );
edgeFace = edgeAttachments( TR, E );
edgeFace = cell2mat( cellfun( resizeCell, edgeFace, ...
    'UniformOutput', false ) );

% #Ex1 boolean list indicating which edges lie on the mesh boundary
bdyEdges = diff(edgeFace) == 0;

% #Fx3 face-edge correspondence tool. Given a list of scalar edge
% quantities, 'EQ', the output of 'EQ(feIDx(f,i))' is that quantity
% corresponding to the edge opposite the ith vertex in face f
e1IDx = sort( [ F(:,3), F(:,2) ], 2 );
e2IDx = sort( [ F(:,1), F(:,3) ], 2 );
e3IDx = sort( [ F(:,2), F(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, sort(E, 2), 'rows' );
[~, e2IDx] = ismember( e2IDx, sort(E, 2), 'rows' );
[~, e3IDx] = ismember( e3IDx, sort(E, 2), 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

%--------------------------------------------------------------------------
% Mesh Geometry Calculations
%--------------------------------------------------------------------------

% Calculate primal mesh edge lengths --------------------------------------
Eij = V(E(:,2), :) - V(E(:,1), :); % Directed primal mesh edge vectors
L_E = sqrt( sum( Eij.^2, 2 ) ); % Defined on edges
L_F = L_E(feIDx); % Defined on faces

% Calculate primal mesh internal angles -----------------------------------

% Some convenience variables to vectorize the cosine law calculation
Gi = L_F; Gj = circshift(L_F, [0 -1]); Gk = circshift(L_F, [0 -2]);

% The internal angles
cosAng = ( Gj.^2 + Gk.^2 - Gi.^2 ) ./ ( 2 .* Gj .* Gk );
sinAng = sin(acos(cosAng));
cotAng = cosAng ./ sinAng;

% Calculate primal mesh face unit normals/areas ---------------------------
e12 = V(F(:,2), :) - V(F(:,1), :);
e13 = V(F(:,3), :) - V(F(:,1), :);

fN = cross( e12, e13, 2 ); 
fA = sqrt( sum( fN.^2, 2 ) ) ./ 2;
fN = fN ./ ( 2 .* fA );

%==========================================================================
% Construct Flat/Sharp Operators
%==========================================================================

% Construct flatPP --------------------------------------------------------
% The 'flatPP' operator maps primal discrete vector fields (tangent vectors
% living on primal vertices) to primal 1-forms living on primal edges. In
% order to work as a single linear operator vector fields should be
% represented as a single (3#V)x1 column vector: [Vx; Vy; Vz]
%
% The 'flatPP' operator can be constructed as the composition of two
% separate operators. The first operator constructs a tangent vector on
% each primal edge by averaging the primal vectors living on the vertices
% defining that edge. The operator is of size (3#E)x(3#V) with 2 non-zero
% elements per row

I = repmat( (1:(3*numE)).', 2, 1 );

J = [ E(:,1); E(:,1)+numV; E(:,1)+(2*numV) ];
J = [ J; E(:,2); E(:,2)+numV; E(:,2)+(2*numV) ];

Q = ones(size(I)) ./ 2;

PP1 = sparse( I, J, Q, 3*numE, 3*numV );

% The second operator calculates the dot product between the tangent
% vectors living on primal edges and the directed edge vector itself,
% producing a single scalar per primal edge. The size of the operator is
% (#E)x(3#E) with 3 non-zero elements per row

I = repmat( (1:numE).', 3, 1 );

J = (1:numE).';
J = [ J; J+numE; J+(2*numE) ];

Q = [ Eij(:,1); Eij(:,2); Eij(:,3) ];

PP2 = sparse(I, J, Q, numE, 3*numE );

this.flatPP = PP2 * PP1;

% Construct flatDP --------------------------------------------------------
% The 'flatDP' operator maps dual discrete vector fields (tangent
% vectors living on dual vertices) to primal 1-forms living on primal
% edges. In order to work as a single linear operator vector fields should
% be represented as a single (3#F)x1 column vector: [ Vx; Vy; Vz ] 
%
% The 'flatDP' operator can be constructed as the composition of two
% separate operators.  The first operator maps calculates the dot product
% of each dual vector with the directed edge vectors defining the boundary
% of the face that contains it.  The operator is of size (3#F)x(3#F) with 3
% non-zero elements per row.  The action of this operator on the tangent
% vector field produces a column vector that if reshaped into the size of
% the face connectivity list would hold for each face the dot product of
% the tangent vector on that face with edge opposite the corresponding
% vertex

I = repmat( (1:(3*numF)).', 3, 1 );

J = repmat( (1:numF).', 3, 1 );
J = [ J; J+numF; J+(2*numF) ];

Q = [ Eij(feIDx(:), 1); Eij(feIDx(:), 2); Eij(feIDx(:), 3) ];

DP1 = sparse( I, J, Q, 3*numF, 3*numF );

% The second operator maps the shadow of the vector fields on the directd
% primal edges into scalar weights defining the primal 1-form field.  It
% does so by constructing a weighted average of the contribution from the
% faces shared by each edge.  The size of the operator is #Ex(3#F)

I = feIDx;

J = (1:numF).';
J = [ J, (J+numF), (J+(2*numF)) ];

% dualLengthInt(f, i) is the length of the intersection of face f with the
% dual edge to the primal edge opposite vertex i in face f.
dualLengthInt = L_F .* cotAng ./ 2;

% dualLength(f, i) is the total length of the dual edge to the primal edge
% opposite vertex i in face f.
dualLength = full( sparse( feIDx, 1, dualLengthInt, numE, 1 ) );
dualLength = dualLength(feIDx);

% The weights are simply the ratio of the intersection of the dual edges in
% the corresponding face to the total length of the dual edge
Q = dualLengthInt ./ dualLength;

DP2 = sparse( I, J, Q, numE, 3*numF );

this.flatDP = DP2 * DP1;

% Construct flatDD --------------------------------------------------------
% The 'flatDD' operator maps dual discrete vector fields (tangent vectors
% living on dual vertices) to dual 1-forms living on dual edges. In order
% to work as a single linear operator vector fields should be represented
% as a single (3#F)x1 column vector: [ Vx; Vy; Vz ]
%
% A single dual edge intersects two triangulation faces (or only one if
% that edge lies on the mesh boundary). The dual 1-form is taken to be the
% sum of the projections of the dual vector fields on each of those faces
% attached to the dual edge with the assocated vector corresponding to the
% intersection of the dual edge with that face.

% Directed edge vectors defined by the orientations of FACES
Ei = V(F(:,3), :) - V(F(:,2), :); % Edge opposite vertex i in face f
Ej = V(F(:,1), :) - V(F(:,3), :); % Edge opposite vertex j in face f
Ek = V(F(:,2), :) - V(F(:,1), :); % Edge opposite vertex k in face f

% The dual edge vectors defined by the orientations of FACES

% Dual edge opposite vertex i in face f
dualEi = cross(fN, Ei ./ repmat(sqrt(sum(Ei.^2, 2)), 1, 3), 2); 
dualEi = repmat(dualLengthInt(:,1), 1, 3) .* dualEi ;

% Dual edge opposite vertex j in face f
dualEj = cross(fN, Ej ./ repmat(sqrt(sum(Ej.^2, 2)), 1, 3), 2); 
dualEj = repmat(dualLengthInt(:,2), 1, 3) .* dualEj ;

% Dual edge opposite vertex k in face f
dualEk = cross(fN, Ek ./ repmat(sqrt(sum(Ek.^2, 2)), 1, 3), 2); 
dualEk = repmat(dualLengthInt(:,3), 1, 3) .* dualEk ;

% We flip the sign of the dual edge where the orientation of the
% face does not match the intrinsic orientation of an edge
si = 1 - 2 .* any( F(:, [2 3]) - E(feIDx(:,1), :), 2 );
sj = 1 - 2 .* any( F(:, [3 1]) - E(feIDx(:,2), :), 2 );
sk = 1 - 2 .* any( F(:, [1 2]) - E(feIDx(:,3), :), 2 );

dualEi = dualEi .* repmat( si, 1, 3 );
dualEj = dualEj .* repmat( sj, 1, 3 );
dualEk = dualEk .* repmat( sk, 1, 3 );

I = repmat(feIDx, 3, 1);

J = (1:numF).';
J = [ J; (J+numF); (J+2*numF) ];
J = repmat(J, 1, 3);

Q = [ dualEi(:), dualEj(:), dualEk(:) ];

this.flatDD = sparse(I, J, Q, numE, 3*numF);

% Construct sharpPD -------------------------------------------------------
% The 'sharpPD' operator maps primal 1-forms living on primal edges into
% dual discrete vector fields (tangent vectors living on dual vertices).
% The operator has size (3#F)x#E.  The output is a single (3#F)x1 column
% vector: [ Vx; Vy; Vz ]
%
% As best I can tell the choice of sharp operator seems to be a little
% ad-hoc.  Here, the sharp operator is basically just an evaluation of a
% continuous 1-form interpolant at each facet centroid.  This choice was
% made so that the gradient operator on primal 0-forms ( sharpDP * d0 )
% matches exactly the classical FEM gradient operator
%
% The interpolation of the 1-forms is performed using the 'Whitney forms'.
%
% The Whitney 0-form, Bi, is simply the familiar 'hat function' defined on
% vertex i, i.e. the unique function that equals one at vertex i, zero at
% all other vertices, and is affine over each 2-simplex
%
% Whitney 1-forms, Wk, interpolate values stored on edges.  For an edge
% oriented away form vertex i towards vertex j: Wk = Bi DBj - Bj DBi, where
% DBi is the gradient of the Whitney 0-form Bi. At the triangle centroid,
% Bi = Bj = Bk = 1/3

% The gradients of the Whitney 0-forms
DBi = cross( fN, Ei, 2 ) ./ ( 2 .* fA );
DBj = cross( fN, Ej, 2 ) ./ ( 2 .* fA );
DBk = cross( fN, Ek, 2 ) ./ ( 2 .* fA );

% The Whitney 1-forms evaluated at face barycenters
Wi = ( DBk - DBj ) ./ 3;
Wj = ( DBi - DBk ) ./ 3;
Wk = ( DBj - DBi ) ./ 3;

% We flip the sign of the Whitney 1-form where the orientation of the
% face does not match the intrinsic orientation of an edge
% si = 1 - 2 .* any( F(:, [2 3]) - E(feIDx(:,1), :), 2 );
% sj = 1 - 2 .* any( F(:, [3 1]) - E(feIDx(:,2), :), 2 );
% sk = 1 - 2 .* any( F(:, [1 2]) - E(feIDx(:,3), :), 2 );

Wi = Wi .* repmat( si, 1, 3 );
Wj = Wj .* repmat( sj, 1, 3 );
Wk = Wk .* repmat( sk, 1, 3 );

I = repmat( (1:(3*numF)).', 3, 1 );

J = repmat( feIDx, 3, 1 );
J = J(:);

Q = [ Wi(:); Wj(:); Wk(:) ];

this.sharpPD = sparse( I, J, Q, 3*numF, numE );

% Construct sharpDD -------------------------------------------------------
% The 'sharpDD' operator maps dual 1-forms living on dual edges into dual
% discrete vector fields (tangent vectors living on dual vertices). The
% operator has size (3#F)x#E. The output is a single (3#F)x1 column vector:
% [Vx; Vy; Vz]
%
% Following Mohamed et al. (2016), the interpolation of dual 1-forms onto
% 2-simplexes is accomplished by interpolating the values of the associated
% primal 1-form using the rotated Whitney 1-forms, i.e. the
% usual Whitney 1-forms rotated CCW 90deg in the plane of the associated
% triangular face.
%
% The 'sharpDD' operator can be constructed as the composition of two
% separate operators. The first operator is simply the Hodge dual
% transforming dual 1-forms into primal 1-forms. This operator is
% constructed separately in the section 'Construct Hodge Dual Operators'.
% The second operator applies the rotated Whitney 1-forms to the associated
% primal 1-forms. Here we simply construct the second operator.

% The signed, rotated Whitney 1-forms
rotWi = cross( fN, Wi, 2 );
rotWj = cross( fN, Wj, 2 );
rotWk = cross( fN, Wk, 2 );

Q = [ rotWi(:); rotWj(:); rotWk(:) ];

incomplete_sharpDD = sparse( I, J, Q, 3*numF, numE );

%==========================================================================
% Construct Exterior Derivative Operators
%==========================================================================

% Construct d0 ------------------------------------------------------------
% The 'd0' operator acts on 0-forms (defined on primal mesh vertices) and
% generates 1-forms (defined on primal mesh edges). #Ex#V sparse matrix

I = repmat( (1:numE).', 2, 1 );
J = E(:);
Q = [ -ones(numE,1); ones(numE,1) ];

this.d0 = sparse( I, J, Q, numE, numV );

% Construct d1 ------------------------------------------------------------
% The 'd1' operator acts on 1-forms (defined on primal mesh edges) and
% generates 2-forms (defined on primal mesh faces). #Fx#E sparse matrix

I = repmat( (1:numF).', 3, 1 );
J = feIDx(:);

% The value of the operator entry is positive where the orientation of the
% edge matches the orientation of the face and negative where the
% orientation of the edge is opposite to the orientation of the face
Q = [ si; sj; sk ];

% Q1 = 1 - 2 .* any( F(:, [2 3]) - E(feIDx(:,1), :), 2 );
% Q2 = 1 - 2 .* any( F(:, [3 1]) - E(feIDx(:,2), :), 2 );
% Q3 = 1 - 2 .* any( F(:, [1 2]) - E(feIDx(:,3), :), 2 );
% 
% Q = [ Q1; Q2; Q3 ];

this.d1 = sparse( I, J, Q, numF, numE );

% Construct dd0 -----------------------------------------------------------
% The 'dd0' operator maps dual 0-forms (defined on dual mesh vertices) and
% to dual 1-forms (defined on dual mesh edges)

this.dd0 = this.d1.';

% Construct dd1 -----------------------------------------------------------
% The 'dd1' operator maps dual 1-forms (defined on dual mesh edges) to dual
% 2-forms (defined on dual mesh faces)

% NOTE: The negative sign is necessary to match a consistent CCW
% orientation of dual 2-cells
this.dd1 = -this.d0.';

%==========================================================================
% Construct Hodge Dual Operators
%==========================================================================
% Here we choose to use the simple diagonal Hodge star operators using the
% circumcentric dual.  Care should be taken when using non-Delaunay meshes
% as the circumcenter of a mesh facet may not lie within its interior.

% Construct hd0 -----------------------------------------------------------
% The 'hd0' operator acts on primal 0-forms (defined on primal mesh
% vertices) and transforms them into dual 2-forms (defined on dual 2-cells)

I = (1:numV).';
J = (1:numV).';

% The weights of the diagonal elements of the 'hd0' operator are simply
% the total area of the dual 2-cells.  These are calculated by determining
% the contribution to this cell from each face adjacent to a vertex and
% then summing these contributions
AV_F = L_F.^2 .* cotAng;
AV_F = ( repmat( sum(AV_F, 2), 1, 3 ) - AV_F ) ./ 8;

Q = sparse( F, 1, AV_F, numV, 1 );

this.hd0 = sparse( I, J, Q, numV, numV );

% Construct hdd2 ----------------------------------------------------------
% The 'hdd2' operator acts on dual 2-forms (defined on dual 2-cells) and
% transforms them into primal 0-forms (defined on primal mesh vertices). In
% the circumcentric dual formulation, it is simply the inverse of the 'hd0'
% operator

this.hdd2 = sparse( I, J, 1./Q, numV, numV );

% Construct hd1 -----------------------------------------------------------
% The 'hd1' operator acts on primal 1-forms (defined on primal mesh edges)
% and tranforms them into dual 1-forms (defined on dual mesh edges)

I = (1:numE).';
J = (1:numE).';

% The weights of the diagonal elements of the 'hd1' operator are simply the
% ratio of the lengths of the dual edge to the corresponding primal edge
Q = sparse( feIDx, 1, cotAng ./ 2, numE, 1 );

this.hd1 = sparse( I, J, Q, numE, numE );

% Assemble the full 'sharpDD' operator
this.sharpDD = incomplete_sharpDD * sparse( I, J, 1./Q, numE, numE );

% Construct hdd1 ----------------------------------------------------------
% The 'hdd1' operator acts on dual 1-forms (defined on dual mesh edges) and
% transforms them into primal 1-forms (defined on primal mesh edges). In
% the circumcentric dual formulation, it is MINUS the inverse of the 'hd1'
% operator (to account for the minus sign accumulated in the operation
% **v = -v for a 1-form v)

this.hdd1 = sparse( I, J, -1./Q, numE, numE );

% Construct hd2 -----------------------------------------------------------
% The 'hd2' operator acts on primal 2-forms (defined on primal mesh faces)
% and transforms them into dual 0-forms (defined on dual mesh vertices)

I = (1:numF).';
J = (1:numF).';
Q = 1 ./ fA;

this.hd2 = sparse( I, J, Q, numF, numF );

% Construct hdd0 ----------------------------------------------------------
% The 'hdd0' operator acts on dual 0-forms (defined on dual mesh vertices)
% and transforms them into primal 2-forms (defined on primal mesh faces).
% In the circumcentric dual formulation, it is simply the inverse of the
% 'hd2' operator

this.hdd0 = sparse( I, J, 1./Q, numF, numF );

end

