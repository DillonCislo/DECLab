classdef DiscreteExteriorCalculus < handle
    %DISCRETEEXTERIORCALCULUS A class used to implement the operators and
    %methods of the Discrete Exterior Calculus on curved/flat mesh
    %triangulations.  See 'Discrete Exterior Calculus' by Matheiu Desbrun
    %et al. (2005) for more details.
    %
    % By Dillon Cislo 01/29/2020
    
    %======================================================================
    %======================================================================
    %                       PROPERTIES
    %======================================================================
    %======================================================================
    
    properties (SetAccess = protected)
        
        % #Fx3 face connectivity list of the input triangulation
        F
        
        % #VxD vertex coordinate list of the input triangulation
        V
        
        % #Ex1 edge connectivity list of the input triangulation
        E
        
        % #Ex#V exterior derivative operator mapping primal 0-forms to
        % primal 1-forms
        d0
        
        % #Fx#E exterior derivative operator mapping primal 1-forms to
        % primal 2-forms
        d1
        
        % #Ex#F exterior derivative operator mapping dual 0-forms to dual
        % 1-forms
        dd0
        
        % #Vx#E exterior derivative operator mapping dual 1-forms to dual
        % 2-forms
        dd1
        
        % #Vx#V Hodge dual operator mapping primal 0-forms to dual 2-forms
        hd0
        
        % #Vx#V Hodge dual operator mapping dual 2-forms to primal 0-forms
        hdd2
        
        % #Ex#E Hodge dual operator mapping primal 1-forms to dual 1-forms
        hd1
        
        % #Ex#E Hodge dual operator mapping dual 1-forms to primal 1-forms
        hdd1
        
        % #Fx#F  Hodge dual operator mapping primal 2-forms to dual 0-forms
        hd2
        
        % #Fx#F Hodge dual operator mapping dual 0-forms to primal 2-forms
        hdd0
        
        % #Ex(3#V) Flat operator mapping primal vector fields to primal
        % 1-forms
        flatPP
        
        % #Ex(3#F) Flat operator mapping dual vector fields to primal
        % 1-forms
        flatDP
        
        % #Ex(3#F) Flat operator mapping dual vector fields to dual 1-forms
        flatDD
        
        % (3#F)x#E Sharp operator mapping primal 1-forms to dual vector
        % fields
        sharpPD
        
        % (3#F)x#E Sharp operator mapping dual 1-forms to dual vector
        % fields
        sharpDD
        
    end
    
    %======================================================================
    %======================================================================
    %                       METHODS
    %======================================================================
    %======================================================================
    methods (Access = public)
        
        function this = DiscreteExteriorCalculus(F, V)
            %DISCRETEEXTERIORCALCULUS Default constructor
            %
            %   INPUT PARAMETERS:
            %
            %       - F:    #Fx3 face connectivity list
            %
            %       - V:    #VxD vertex coordinate list
            %
            %   OUTPUT PARAMETERS:
            %
            %       - this: An instance of this class
            
            % Validate Inputs ---------------------------------------------
            validateattributes( V, {'numeric'}, ...
                {'2d', 'finite', 'nonnan', 'real'} );
            validateattributes( F, {'numeric'}, ...
                {'2d', 'ncols', 3, 'finite', 'nonnan', ...
                'integer', 'real', 'positive'} );
            
            assert( all(ismember(F(:), (1:size(V,1)).')), ...
                'Face connectivity list contains unreferenced vertices' );
            
            this.F = F;
            this.V = V;
            this.E = sort(edges( triangulation(F, V) ), 2);
            
            % Construct DEC Operators -------------------------------------
            this.constructDECOperators();
            
        end
        
        function UFlat = primalVectorToPrimal1Form(this, U)
            %PRIMALVECTORTOPRIMAL1FORM Transforms the tangent part of a
            %vector field living on triangulation vertices (primal
            %vertices) into a primal 1-form field. An implementation of a
            %'discrete flat' operator. Vertex normals are estimated using
            %built-in MATLAB face area weighting. This method assumes that
            %the face connectivity list is consistently ordered to produce
            %consistent vertex normals
            
            % Input Processing --------------------------------------------
            validateattributes( U, {'numeric'}, ...
                {'2d', 'finite', 'nonnan', 'real', ...
                'ncols', size(this.V, 2), 'nrows', size(this.V, 1)} );
            
            % Project the vector field onto the tangent space
            % of its respective vertex ------------------------------------
            VN = vertexNormal(triangulation(this.F, this.V));
            U = U - repmat(dot(U, VN, 2), 1, 3) .* VN;
            
            % Apply flat operator to construct primal 1-form --------------
            UFlat = this.flatPP * U(:);
            
        end
        
        function UFlat = dualVectorToPrimal1Form(this, U)
            %DUALVECTORTOPRIMAL1FORM Transforms the tangent part of a
            %vector field living on triangulation faces (dual vertices)
            %into a primal 1-form field.  An implementation of a 'discrete
            %flat' operator
            %
            %   INPUT PARAMETERS:
            %
            %       - U:        #FxD vector field
            %
            %   OUTPUT PARAMETERS:
            %
            %       - UFlat:    #Ex1 primal 1-form
            
            % Input Processing --------------------------------------------
            validateattributes( U, {'numeric'}, ...
                {'2d', 'finite', 'nonnan', 'real', ...
                'ncols', size(this.V, 2), 'nrows', size(this.F, 1)} );
            
            % Project the vector field onto the tangent space of its
            % respective face ---------------------------------------------
            FN = faceNormal(triangulation(this.F, this.V));
            U = U - repmat(dot(U, FN, 2), 1, 3) .* FN;
            
            % Apply flat operator to construct primal 1-form --------------
            UFlat = this.flatDP * U(:);
            
        end
        
        function UFlat = dualVectorToDual1Form(this, U)
            %DUALVECTORTODUAL1FORM Transforms the tangent part of a vector
            %field living on triangulation faces (dual vertices) into a
            %dual 1-form field. An implementation of a 'discrete flat'
            %operator
            %
            %   INPUT PARAMETERS:
            %
            %       - U:        #FxD vector field
            %
            %   OUTPUT PARAMETERS:
            %
            %       - UFlat:    #Ex1 dual 1-form
            
            % Input Processing --------------------------------------------
            validateattributes( U, {'numeric'}, ...
                {'2d', 'finite', 'nonnan', 'real', ...
                'ncols', size(this.V, 2), 'nrows', size(this.F, 1)} );
            
            % Project the vector field onto the tangent space of its
            % respective face ---------------------------------------------
            FN = faceNormal(triangulation(this.F, this.V));
            U = U - repmat(dot(U, FN, 2), 1, 3) .* FN;
            
            % Apply flat operator to construct dual 1-form ----------------
            UFlat = this.flatDD * U(:);
            
        end
            
        
        function USharp = primal1FormToDualVector(this, U)
            %PRIMAL1FORMTODUALVECTOR Transforms a primal 1-form field
            %living on mesh edges into a tangent dual vector field living
            %on triangulation faces (dual vertices)
            %
            %   INPUT PARAMETERS:
            %
            %       - U:        #Ex1 primal 1-form
            %
            %   OUTPUT PARAMETERS:
            %
            %       - USharp:   #FxD dual vector field
            
            % Input Processing --------------------------------------------
            validateattributes( U, {'numeric'}, ...
                {'2d', 'vector', 'numel', size(this.E, 1), 'finite', ...
                'real', 'nonnan'} );
            
            % Ensure that U is a column vector
            if (size(U, 2) ~= 1), U = U.'; end
            
            % Apply sharp operator to construct dual vector field ---------
            USharp = this.sharpPD * U;
            
            % Re-shape to appropriate output dimensions
            USharp = reshape( USharp, size(this.F, 1), size(this.V, 2) );
            
            % Project the vector field onto the tangent space of its
            % respective face ---------------------------------------------
            FN = faceNormal(triangulation(this.F, this.V));
            USharp = USharp - repmat(dot(USharp, FN, 2), 1, 3) .* FN;
            
        end
        
        function USharp = dual1FormToDualVector(this, U)
            %DUAL1FORMTODUALVECTOR Transforms a dual 1-form field living on
            %mesh (dual) edges into a tangent dual vector field living on
            %triangulation faces (dual vertices)
            %
            %   INPUT PARAMETERS:
            %
            %       - U:        #Ex1 dual 1-form
            %
            %   OUTPUT PARAMETERS:
            %
            %       - USharp:   #FxD dual vector field
            
            % Input Processing --------------------------------------------
            validateattributes( U, {'numeric'}, ...
                {'2d', 'vector', 'numel', size(this.E, 1), 'finite', ...
                'real', 'nonnan'} );
            
            % Ensure that U is a column vector
            if (size(U, 2) ~= 1), U = U.'; end
            
            % Apply sharp operator to construct dual vector field ---------
            USharp = this.sharpDD * U;
            
            % Re-shape to appropriate output dimensions
            USharp = reshape( USharp, size(this.F, 1), size(this.V, 2) );
            
            % Project the vector field onto the tangent space of its
            % respective face ---------------------------------------------
            FN = faceNormal(triangulation(this.F, this.V));
            USharp = USharp - repmat(dot(USharp, FN, 2), 1, 3) .* FN;
            
        end
        
        function gradS = gradient(this, S)
            %GRADIENT Calculates the discrete gradient of a primal 0-form
            %(scalar function defined on mesh vertices)
            %
            %   INPUT PARAMETERS:
            %
            %       - S:        #Vx1 primal 0-form
            %
            %   OUTPUT PARAMETERS:
            %
            %       - gradS:    #Fx3 dual vector field
            
            % Input Processing --------------------------------------------
            validateattributes( S, {'numeric'}, ...
                {'2d', 'vector', 'numel', size(this.V, 1), 'finite', ...
                'real', 'nonnan'} );
            
            % Ensure that S is a column vector
            if (size(S, 2) ~= 1), S = S.'; end
            
            % Calculate gradient ------------------------------------------
            gradS = this.primal1FormToDualVector(this.d0 * S);
            
        end
        
        function divU = divergence(this, U, route)
            %DIVERGENCE Calculates the discrete divergence of a dual
            %tangent vector field or something like a divergence of a
            %primal/dual 1-form
            %
            %   INPUT PARAMETERS:
            %
            %       - U:        #FxD dual vector field OR
            %                   #VxD primal vector field OR
            %                   #Ex1 primal 1-form OR
            %                   #Ex1 dual 1-form
            %
            %       - route:    The type of route to take when calculating
            %                   the divergence
            %                   'primal' (default) maps a primal 1-form to
            %                   a primal 0-form
            %                   'dual' maps a dual 1-form to a dual 0-form
            %
            %   OUTPUT PARAMETERS
            %
            %       - divU:     #Vx1 primal 0-form OR
            %                   #Fx1 dual 0-form
            
            % Input Processing --------------------------------------------
            if (nargin < 3), route = 'primal'; end
            
            validateattributes( U, {'numeric'}, ...
                {'2d', 'finite', 'real', 'nonnan'} );
            
            if ~(strcmpi(route, 'primal') || strcmpi(route, 'dual'))
                error('Invalid primal/dual calculation route');
            end
            
            if ~isvector(U)
                
                % If the input field is a tangent primal/dual vector field
                % we must convert it into a primal 1-form
                assert(size(U, 2) == size(this.V, 2), ...
                    'Input vector field has incorrect dimensions!');
                
                if (size(U, 1) == size(this.F, 1))
                    
                    if strcmpi(route, 'primal')
                        U = this.dualVectorToPrimal1Form(U);
                    else
                        U = this.dualVectorToDual1Form(U);
                    end
                    
                elseif (size(U, 1) == size(this.V, 1))
                    
                    if strcmpi(route, 'primal')
                        U = this.primalVectorToPrimal1Form(U);
                    else
                        error(['A primal->dual flat operator has not ' ...
                            'yet been implemented']);
                    end
                    
                else
                    
                    error('Input vector field is improperly sized!');
                    
                end
                
            else
                
                assert( numel(U) == size(this.E, 1), ...
                    'Input primal 1-form is improperly sized!');
                
                if (size(U,2) ~= 1), U = U.'; end
                
            end
            
            % Calculate divergence ----------------------------------------
            
            if strcmpi(route, 'primal')
                % divU = inv(this.hd0) * this.dd1 * this.hd1 * U;
                divU = this.hdd2 * this.dd1 * this.hd1 * U;
            else
                % divU = this.hd2 * this.d1 * (-inv(this.hd1)) * U;
                divU = this.hd2 * this.d1 * this.hdd1 * U;
            end
            
        end
        
        function curlU = curl(this, U, route)
            %CURL Calculates something like a discrete curl for a
            %primal/dual tangent vector field or a primal/dual 1-form.
            %
            %   INPUT PARAMETERS:
            %
            %       - U:        #FxD dual vector field OR
            %                   #VxD primal vector field OR
            %                   #Ex1 primal 1-form OR
            %                   #Ex1 dual 1-form
            %
            %       - route:    The type of route to take when calculating
            %                   the divergence
            %                   'primal' (default) maps a primal 1-form to
            %                   a dual 0-form
            %                   'dual' maps a dual 1-form to a primal
            %                   0-form
            %
            %   OUTPUT PARAMETERS
            %
            %       - curlU:    #Fx1 dual 0-form OR
            %                   #Vx1 primal 0-form
          
            % Input Processing --------------------------------------------
            if (nargin < 3), route = 'primal'; end
            
            validateattributes( U, {'numeric'}, ...
                {'2d', 'finite', 'real', 'nonnan'} );
            
            if ~(strcmpi(route, 'primal') || strcmpi(route, 'dual'))
                error('Invalid primal/dual calculation route');
            end
            
            if ~isvector(U)
                
                % If the input field is a tangent primal/dual vector field
                % we must convert it into a primal 1-form
                assert(size(U, 2) == size(this.V, 2), ...
                    'Input vector field has incorrect dimensions!');
                
                if (size(U, 1) == size(this.F, 1))
                    
                    if strcmpi(route, 'primal')
                        U = this.dualVectorToPrimal1Form(U);
                    else
                        U = this.dualVectorToDual1Form(U);
                    end
                    
                elseif (size(U, 1) == size(this.V, 1))
                    
                    if strcmpi(route, 'primal')
                        U = this.primalVectorToPrimal1Form(U);
                    else
                        error(['A primal->dual flat operator has not ' ...
                            'yet been implemented']);
                    end
                    
                else
                    
                    error('Input vector field is improperly sized!');
                    
                end
                
            else
                
                assert( numel(U) == size(this.E, 1), ...
                    'Input primal 1-form is improperly sized!');
                
                if (size(U,2) ~= 1), U = U.'; end
                
            end
            
            % Calculate curl ----------------------------------------------
            
            if strcmpi(route, 'primal')
                curlU = this.hd2 * this.d1 * U;
            else
                % curlU = inv(this.hd0) * this.dd1 * U;
                curlU = this.hdd2 * this.dd1 * U;
            end
            
        end
        
        function [ divU, rotU, harmU, scalarP, vectorP ] = ...
                helmholtzHodgeDecomposition(this, U, r)
            %HELMHOLTZHODGEDECOMPOSITION Calculate the Helmholtz-Hodge
            %decomposition of a tangent dual vector field or a primal
            %1-form.
            %
            %   INPUT PARAMETERS:
            %
            %       - U:        #FxD dual vector field OR
            %                   #Ex1 primal 1-form
            %
            %       - r:        Small number used for Tikhonov
            %                   regularization of the symmetric Laplace
            %                   operators prior to linear solves
            %
            %   OUTPUT PARAMETERS:
            %
            %       NOTE: The size of the first three outputs will match
            %       the size of the input argument
            %
            %       - divU:     The irrotational (curl-free) part of U
            %
            %       - rotU:     The divergence-free part of U
            %       
            %       - harmU:    The harmonic part of U
            %
            %       - scalarP:  #Vx1 scalar potential of U (primal 0-form)
            %
            %       - vectorP:  #Fx1 vector potential of U (primal 2-form)
            
            % Input Processing --------------------------------------------
            validateattributes( U, {'numeric'}, ...
                {'2d', 'finite', 'real', 'nonnan'} );
            
            if (nargin < 3), r = 0; end
            validateattributes( r, {'numeric'}, ...
                {'scalar', 'nonnegative', 'finite', 'real', 'nonnan'});
            
            if ~isvector(U)
                
                % If the input field is a tangent dual vector field we must
                % convert it into a primal 1-form
                assert( (size(U, 1) == size(this.F, 1)) &&...
                    (size(U, 2) == size(this.V, 2)), ...
                    'Input dual vector field is improperly sized!');
                
                UVec = U; % Store a copy of input vector field
                
                U = this.dualVectorToPrimal1Form(U);
                
                isVector = true;
                
            else
                
                assert( numel(U) == size(this.E, 1), ...
                    'Input primal 1-form is improperly sized!');
                
                if (size(U,2) ~= 1), U = U.'; end
                
                isVector = false;
                
            end
            
            % Peform the decomposition ------------------------------------
            
            % Calculate the scalar potential. Use of double negatives here
            % makes the operators symmetric positive semi-definite
            LS = -this.dd1 * this.hd1 * this.d0;
            if (r > 0), LS = LS + r * speye(size(LS)); end
            scalarP = LS \ ( -this.dd1 * this.hd1 * U );
            divU = this.d0 * scalarP;
            
            % Calculate the vector potential. Use of double negatives here
            % makes the operators symmetric positive semi-definite
            % (MINUS SIGNS SHOULD NOW BE ACCOUNTED FOR)
            
            % OLD WAY: 'rotU' was the same, but 'vectorP' was off by a
            % minus sign
            % vectorP = ( this.d1 * inv(this.hd1) * this.dd0 ) \ ...
            %     ( this.d1 * U );
            % vectorP = inv(this.hd2) * vectorP;
            % rotU = inv(this.hd1) * this.dd0 * this.hd2 * vectorP;
            
            % NEW WAY
            LV = -this.d1 * this.hdd1 * this.dd0;
            if (r > 0), LV = LV + r * speye(size(LV)); end
            vectorP = LV \ (-this.d1 * U );
            vectorP = this.hdd0 * vectorP;
            rotU = this.hdd1 * this.dd0 * this.hd2 * vectorP;
            
            % Calculate the harmonic part of U
            if isVector
                
                divU = this.primal1FormToDualVector(divU);
                try
                    rotU = this.primal1FormToDualVector(rotU);
                catch
                    debugMsg(1, 'Cannot compute rotU, setting to NaN \n')
                    rotU = NaN * ones(size(divU)) ;
                end
                harmU = UVec - divU - rotU;
                
            else
                
                harmU = U - divU - rotU;
                
            end
            
        end
        
        function lapU = laplacian(this, U, normalizeAreas)
            %LAPLACIAN Calculate the action of the discrete
            %Laplace-Beltrami operator on an input field
            %
            %   INPUT PARAMETERS:
            %
            %       - U:        #Vx1 primal 0-form OR
            %                   #Vx3 primal tangent vector field OR
            %                   #Fx3 dual tangent vector field OR
            %                   #Ex1 primal 1-form OR
            %                   #Fx1 primal 2-form
            %
            %       - normalizeAreas:   A boolean. If true the Laplacian
            %       operator will be normalized by the area of dual
            %       2-simplices (i.e. scalar Laplacian output is a primal
            %       0-form).  If false, the Laplacian operator will NOT be
            %       normalized by these areas (i.e. scalar Laplacian output
            %       is a dual 2-form).  For reference, the latter
            %       formalism corresponds exactly to the classical FEM
            %       cotangent Laplacian with no area weights.
            %
            %   OUTPUT PARAMETERS:
            %
            %       NOTE: The size of the output will depend on the
            %       the size of the input argument
            %
            %       - lapU:     The Laplacian of the input field
            
            % Input Processing --------------------------------------------
            validateattributes( U, {'numeric'}, ...
                {'2d', 'finite', 'real', 'nonnan'} );
            
            if size(U, 1) == size(this.V, 1)
                
                if isvector(U)
                    
                    inputType = '0form';
                    
                else
                    
                    error(['This method cannot handle primal tangent ' ...
                        'vector inputs yet']);
                    
                end
            
            elseif ~isvector(U)
                
                % If the input field is a tangent dual vector field we must
                % convert it into a primal 1-form
                assert( (size(U, 1) == size(this.F, 1)) &&...
                    (size(U, 2) == size(this.V, 2)), ...
                    'Input dual vector field is improperly sized!');
                
                inputType = 'dualVector';
                
            else
                
                if numel(U) == size(this.E, 1)
                    
                    error(['Direct 1-form handling not yet ' ...
                        'implemented for this method']);
                    inputType = '1form';
                    
                elseif numel(U) == size(this.V, 1)
                    
                    inputType = '0form';
                    
                elseif numel(U) == size(this.F, 1)
                    
                    inputType = '2form';
                    
                else
                    
                    error('Input field is improperly sized!');
                    
                end
                
                if (size(U,2) ~= 1), U = U.'; end
                
            end
            
            % Set default for area normalization
            if (nargin < 3), normalizeAreas = true; end
            
            % Calculate the Laplacian of the input field ------------------
            
            switch inputType
                
                case '0form'
                    
                    if normalizeAreas
                        
                        % lapU = inv(this.hd0) * this.dd1 * ...
                        %     this.hd1 * this.d0 * U;
                        lapU = this.hdd2 * this.dd1 * ...
                            this.hd1 * this.d0 * U;
                    else
                        
                        lapU = this.dd1 * this.hd1 * this.d0 * U;
                        
                    end
                    
                case '1form'
                    
                    error(['Direct 1-form handling not yet ' ...
                        'implemented for this method']);
                    
                    % Primal route (d * d * u works, but * d * d u does not)
                    % lapU = ( this.hdd1 * this.dd0 * this.hd2 * this.d1 + ...
                    %     this.d0 * this.hdd2 * this.dd1 * this.hd1 ) * U;
                    
                    % Dual route (* d * d u works, but d * d * u does not)
                    lapU = ( this.hd1 * this.d0 * this.hdd2 * this.dd1 + ...
                        this.dd0 * this.hd2 * this.d1 * this.hdd1 ) * U;
                    
                case 'dualVector'
                    
                    % d * d * u
                    lapU1 = this.dualVectorToPrimal1Form(U);
                    lapU1 = this.d0 * this.hdd2 * this.dd1 * this.hd1 * lapU1;
                    lapU1 = this.primal1FormToDualVector(lapU1);
                    
                    % * d * d u
                    lapU2 = this.dualVectorToDual1Form(U);
                    lapU2 = this.hd1 * this.d0 * this.hdd2 * this.dd1 * lapU2;
                    lapU2 = this.dual1FormToDualVector(lapU2);
                    
                    lapU = lapU1 + lapU2;
                    
                case '2form'
                    
                    lapU = this.d1 * this.hdd1 * this.dd0 * this.hd2 * U;
                    
            end
            
        end
        
        function Fnew = CCWOrientFaces(this, F, V)
            %CCWORIENTFACES Re-orders the face connectivty list of an input
            %mesh triangulation so that all faces are counter-clockwise
            %ordered in the plane of the face
            %
            %   INPUT PARAMETERS:
            %
            %       - F:        #Fx3 face connectivity list
            %
            %       - V:        #VxD vertex coordinate list
            %
            %   OUTPUT PARAMETERS:
            % 
            %       - Fnew:     #Fx3 re-ordred face connectivity list
            
            % Validate Inputs ---------------------------------------------
            validateattributes( V, {'numeric'}, ...
                {'2d', 'finite', 'nonnan', 'real'} );
            validateattributes( F, {'numeric'}, ...
                {'2d', 'ncols', 3, 'finite', 'nonnan', ...
                'integer', 'real', 'positive'} );
            
            assert( all(ismember(F(:), (1:size(V,1)).')), ...
                'Face connectivity list contains unreferenced vertices' );
            
            % Determine current face orientation --------------------------
            
            % The extrinsic edge vector for the edge opposite vertex k
            Ek = V(F(:,2), :) - V(F(:,1), :);
            Lk = sqrt( sum( Ek.^2, 2 ) ); % The length of the edge
            
            % The extrinsic edge vector for the edge opposite vertex j
            Ej = V(F(:,3), :) - V(F(:,1), :);
            Lj = sqrt( sum( Ej.^2, 2 ) ); % The length of the edge
            
            % The length of the edge vector opposite vertex i
            Li = V(F(:,3), :) - V(F(:,2), :);
            Li = sqrt( sum( Li.^2, 2 ) );
            
            % Functions of the internal angle associate to vertex i
            cosAng = ( Lj.^2 + Lk.^2 - Li.^2 ) ./ ( 2 .* Lj .* Lk );
            sinAng = sin(acos(cosAng));
            
            % The intrinsic edge vector for the edge opposite vertex k
            ek = [ Lk, zeros(size(F,1), 2) ];
            
            % The intrinsic edge vector for the edge opposite vertex j
            ej = Lj .* [ cosAng, sinAng, zeros(size(F,1), 1) ];
            
            % Faces are CCW ordered if the sign of the 3rd component of the
            % cross product of the intrinsic edges are positive
            faceOrder = cross(ek, ej, 2);
            faceOrder = sign( faceOrder(:, 3) );
            
            % Update face orientation -------------------------------------
            
            Fnew = F;
            for f = 1:size(F,1)
                
                if faceOrder(f) < 0
                    Fnew(f,:) = fliplr(Fnew(f,:));
                end
                
            end
            
        end

    end

end

