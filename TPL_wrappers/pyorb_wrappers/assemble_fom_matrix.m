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
        [ A, ~ ] = assembler_poisson( fespace, f, mu,dirichlet_functions, neumann_functions );
        if(isfield(fem_specifics,'final_time'))
            M = assemble_mass(fespace);
        end
    else
        element_list = varargin{1};
        [ A, ~ ] = assembler_poisson( fespace, f, mu,dirichlet_functions, neumann_functions, element_list );
        if(isfield(fem_specifics,'final_time'))
            M = assemble_mass(fespace,1,element_list);
        end
    end
    
    if nargin < 4
        
       [ i, j, val ] = find( A );
       array.A = [ i, j, val ];
       if(isfield(fem_specifics,'final_time'))
           
           [ i, j, val ] = find( M );
           array.M = [ i, j, val ];

       end
       
    elseif nargin == 4
        
        indices_list = varargin{2}; % supposed to be a matrix of size nb_of_indices x 2 
        elements_A = zeros( size( indices_list, 1 ), 1 );
        for ii = 1:size( indices_list , 1 )
            elements_A(ii) = A(indices_list(ii, 1), indices_list(ii, 2));
        end       
        array.A = [indices_list, elements_A];
        
        if(isfield(fem_specifics,'final_time'))
            elements_M = zeros( size( indices_list, 1 ), 1 );
        for ii = 1:size( indices_list , 1 )
            elements_M(ii) = M(indices_list(ii, 1), indices_list(ii, 2));
        end       
        array.M = [indices_list, elements_M];    
        end
        
    end
    

    
end

