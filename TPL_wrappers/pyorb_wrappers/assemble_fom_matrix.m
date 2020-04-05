function [array] = assemble_fom_matrix( param, fem_specifics, varargin )
% Assemble FEM approximation of the matrices A (stiffness) and M (mass). Also, 
% it is possible to evaluate it just at some mesh elements or DOF indices via 
% varagin{1} and varagin{2} respectively.
% input=
%           param: vector of characteristic parameters; if 1, default stiffness matrix is
%           constructed
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
%           varargin: if passed, it contains as first element the list of
%           mesh elements to be considered in computing the matrices and as
%           second element the list of DOF indices at which evaluating the
%           matrices
% output=
%           array: struct containing the mass and stiffness matrices in COO format

    [~, fespace] = set_fem_simulation( fem_specifics, [1,1,0,1]);

    dirichlet_functions = @(x) [0;0;0;0];
    neumann_functions = @(x) [0;0;0;0];

    % forcing term (not employed)
    f = @(x) 0*x(1,:);

    current_model = fem_specifics.model;
    
    if param == 1
        mu = param;
    else
        mu = build_diffusion( param, current_model );
    end
    
    if nargin == 2
        [ A, ~ ] = assembler_poisson( fespace, f, mu, dirichlet_functions, neumann_functions );
        M = assemble_mass(fespace);
    else
        element_list = varargin{1};
        [ A, ~ ] = assembler_poisson( fespace, f, mu, dirichlet_functions, neumann_functions, element_list );
        M = assemble_mass(fespace,1,element_list);
    end
    
    if nargin < 4
        
       [ i, j, val ] = find( A );
       array.A = [ i, j, val ];
           
       [ i, j, val ] = find( M );
       array.M = [ i, j, val ];
       
    elseif nargin == 4
        
        indices_list = varargin{2}; % supposed to be a matrix of size nb_of_indices x 2 
        elements_A = zeros( size( indices_list, 1 ), 1 );
        for ii = 1:size( indices_list , 1 )
            elements_A(ii) = A(indices_list(ii, 1), indices_list(ii, 2));
        end       
        array.A = [indices_list, elements_A];
        elements_M = zeros( size( indices_list, 1 ), 1 );
        for ii = 1:size( indices_list , 1 )
            elements_M(ii) = M(indices_list(ii, 1), indices_list(ii, 2));
        end       
        array.M = [indices_list, elements_M];   
        
    end
    
end

