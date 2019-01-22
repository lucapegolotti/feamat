function [array] = assemble_fom_matrix( param, fem_specifics, varargin )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
% output=
%           array: struct containing the stiffness matrix in COO format

    [~, fespace] = set_fem_simulation( fem_specifics );

    dirichlet_functions = @(x) [0;0;0;0];
    neumann_functions = @(x) [1;0;0;0];

    % forcing term (not employed)
    f = @(x) 0*x(1,:);

    current_model = fem_specifics.model;

    mu = build_diffusion( param, current_model );
    
    if nargin == 2
        [ A, b ] = assembler_poisson( fespace, f, mu,dirichlet_functions, neumann_functions );
    else
        element_list = varargin{1};
        [ A, b ] = assembler_poisson( fespace, f, mu,dirichlet_functions, neumann_functions, element_list );
    end
    
    if nargin < 4
       [ i, j, val ] = find( A );
       array.A = [ i, j, val ];
    else if nargin == 4
        indeces_list = varargin{2}; % supposed to be a matrix of size nb_of_indices x 2 
        elements_A = zeros( size( indeces_list, 1 ), 1 );
        for ii = 1:size( indeces_list , 1 )
            elements_A(ii) = A(indeces_list(ii, 1), indeces_list(ii, 2));
        end
        
        array.A = [indeces_list, elements_A];
    end
    

    
end

