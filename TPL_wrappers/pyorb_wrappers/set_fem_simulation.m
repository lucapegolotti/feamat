function [mesh, fespace] = set_fem_simulation( fem_specifics, varargin )
% Assemble fom affine matrix for elliptic scalar problems
% input=
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           mesh: mesh for the FE problem to be solved in feamat
%           fespace: fespace for the FE problem to be solved in feamat
 
n_elements_x = fem_specifics.number_of_elements;
poly_degree  = fem_specifics.polynomial_degree;

n_elements_y = n_elements_x;

% geometry  definition
bottom_left_corner_x = 0;
bottom_left_corner_y = 0;

L = 1.0;
H = 1.0;

file_name = '';

if nargin > 1
    file_name = varargin{1};
end

mesh = create_mesh( bottom_left_corner_x, ...
                    bottom_left_corner_y, ...
                    L,H,n_elements_x,n_elements_y, ... 
                    file_name );

current_model = fem_specifics.model;
                      
bc_flags = [0 0 1 0];

if strcmp( current_model, 'nonaffine' )
    bc_flags = [1 1 1 1];
end

current_dirichlet = fem_specifics.use_nonhomogeneous_dirichlet;
if strcmp( current_model, 'nonaffine' ) && strcmp( current_dirichlet, 'Y' )
    bc_flags = [1 1 0 1];
end

fespace = create_fespace( mesh, poly_degree, bc_flags );

end