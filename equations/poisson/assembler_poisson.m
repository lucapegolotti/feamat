function [A,b] = assembler_poisson(fespace,fun,mu,dirichlet_functions,neumann_functions, varargin )
% Assemble poisson matrix with boundary conditions
% input=
%           fespace: finite elemnet space
%           fun: anonymous function of the forcing term
%           mu: anonymous function or scalar of the diffusion coefficient
%               If scalar, the code is optimized on structured meshes
%           dirichlet_functions: Dirichlet boundary data
%           neumann_functions: Neumann_boundary data
% output=
%           A: system matrix
%           b: right handside

use_full_element_list = true;

if nargin > 5
    use_full_element_list = false;
    element_list = varargin{1};
end

bc_flags = fespace.bc;

thereisneumann = 1;

if (length(find(bc_flags)) == length(bc_flags))
    thereisneumann = 0;
end

if use_full_element_list
    A = assemble_stiffness( mu, fespace );
else
    A = assemble_stiffness_elementlist( mu, fespace, element_list );
end

b = assemble_rhs(fespace,fun);

if (thereisneumann)
   b = apply_neumann_bc(b,fespace,neumann_functions); 
end

% Apply Dirichlet boundary conditions
[A,b] = apply_dirichlet_bc(A,b,fespace,dirichlet_functions);


