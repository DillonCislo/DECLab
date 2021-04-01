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
        
        % #Vx#E exterior derivative operator mapping primal 0-forms to
        % primal 1-forms
        d0
        
        % #Ex#F exterior derivative operator mapping primal 1-forms to
        % primal 2-forms
        d1
        
        % #Fx#E exterior derivative operator mapping dual 0-forms to dual
        % 1-forms
        dd0
        
        % #Ex#V exterior derivative operator mapping dual 1-forms to dual
        % 2-forms
        dd1
        
        % #Vx#V Hodge dual operator mapping primal 0-forms to dual 2-forms
        hd0
        
        % #Ex#E Hodge dual operator mapping primal 1-forms to dual 1-forms
        hd1
        
        % #Fx#F  Hodge dual operator mapping primal 2-forms to dual 0-forms
        hd2
        
        % #Ex(3#F) Flat operator mapping dual vector fields to primal
        % 1-forms
        flatDP
        
        % (3#F)x#E Sharp operator mapping primal 1-forms to dual vector
        % fields
        sharpPD
        
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
            this.E = edges( triangulation(F, V) );
            
            % Construct DEC Operators -------------------------------------
            this.constructDECOperators();
            
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
            
            % Apply flat operator to construct primal 1-form --------------
            UFlat = this.flatDP * U(:);
            
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
            gradS = this.sharpPD * this.d0 * S;
            
            % Re-shape to appropriate output dimensions
            gradS = reshape( gradS, size(this.F, 1), size(this.V, 2) );
            
        end
        
        function divU = divergence(this, U)
            %DIVERGENCE Calculates the discrete divergence of a dual
            %tangent vector field or something like a divergence of a
            %primal 1-form
            %
            %   INPUT PARAMETERS:
            %
            %       - U:    #FxD dual vector field OR
            %               #Ex1 primal 1-form
            %
            %   OUTPUT PARAMETERS
            %
            %       - divU: #Vx1 primal 0-form
            
            % Input Processing --------------------------------------------
            validateattributes( U, {'numeric'}, ...
                {'2d', 'finite', 'real', 'nonnan'} );
            
            if ~isvector(U)
                
                % If the input field is a tangent dual vector field we must
                % convert it into a primal 1-form
                assert( (size(U, 1) == size(this.F, 1)) &&...
                    (size(U, 2) == size(this.V, 2)), ...
                    'Input dual vector field is improperly sized!');
                
                U = this.dualVectorToPrimal1Form(U);
                
            else
                
                assert( numel(U) == size(this.E, 1), ...
                    'Input primal 1-form is improperly sized!');
                
                if (size(U,2) ~= 1), U = U.'; end
                
            end
            
            % Calculate divergence ----------------------------------------
            divU = inv(this.hd0) * this.dd1 * this.hd1 * U;
            
        end
        
        function curlU = curl(this, U)
            %CURL Calculates something like a discrete curl for a dual
            %tangent vector field or a primal 1-form.
            %
            %   INPUT PARAMETERS:
            %
            %       - U:    #FxD dual vector field OR
            %               #Ex1 primal 1-form
            %
            %   OUTPUT PARAMETERS
            %
            %       - curlU: #Fx1 dual 0-form
            
            % Input Processing --------------------------------------------
            validateattributes( U, {'numeric'}, ...
                {'2d', 'finite', 'real', 'nonnan'} );
            
            if ~isvector(U)
                
                % If the input field is a tangent dual vector field we must
                % convert it into a primal 1-form
                assert( (size(U, 1) == size(this.F, 1)) &&...
                    (size(U, 2) == size(this.V, 2)), ...
                    'Input dual vector field is improperly sized!');
                
                U = this.dualVectorToPrimal1Form(U);
                
            else
                
                assert( numel(U) == size(this.E, 1), ...
                    'Input primal 1-form is improperly sized!');
                
                if (size(U,2) ~= 1), U = U.'; end
                
            end
            
            % Calculate curl ----------------------------------------------
            curlU = this.hd2 * this.d1 * U;
            
        end
        
        function [ divU, rotU, harmU, scalarP, vectorP ] = ...
                helmholtzHodgeDecomposition(this, U)
            %HELMHOLTZHODGEDECOMPOSITION Calculate the Helmholtz-Hodge
            %decomposition of a tangent dual vector field or a primal
            %1-form.
            %
            %   INPUT PARAMETERS:
            %
            %       - U:        #FxD dual vector field OR
            %                   #Ex1 primal 1-form
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
            % WHAT THE FUCK IS WITH THESE MINUS SIGNS
            
            % Calculate the scalar potential
            scalarP = ( this.dd1 * this.hd1 * this.d0 ) \ ...
                ( this.dd1 * this.hd1 * U );
            divU = this.d0 * scalarP;
            
            % Calculate the vector potential (SIGN CHANGE CHECK)
            vectorP = ( this.d1 * inv(this.hd1) * this.dd0 ) \ ...
                ( this.d1 * U );
            vectorP = inv(this.hd2) * vectorP;
            rotU = inv(this.hd1) * this.dd0 * this.hd2 * vectorP;
            
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
            %                   #Fx3 dual tangent vector field OR
            %                   #Ex1 primal 1-form
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
                
                inputType = '0form';
            
            elseif ~isvector(U)
                
                % If the input field is a tangent dual vector field we must
                % convert it into a primal 1-form
                assert( (size(U, 1) == size(this.F, 1)) &&...
                    (size(U, 2) == size(this.V, 2)), ...
                    'Input dual vector field is improperly sized!');
                
                U = this.dualVectorToPrimal1Form(U);
                
                inputType = 'dualVector';
                
            else
                
                if numel(U) == size(this.E, 1)
                    
                    inputType = '1form';
                    
                elseif numel(U) == size(this.V, 1)
                    
                    inputType = '0form';
                    
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
                        lapU = inv(this.hd0) * this.dd1 * ...
                            this.hd1 * this.d0 * U;
                    else
                        lapU = this.dd1 * this.hd1 * this.d0 * U;
                    end
                    
                case '1form'
                    
                    lapU = [];
                    
                case 'dualVector'
                    
                    lapU = [];
                    
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

