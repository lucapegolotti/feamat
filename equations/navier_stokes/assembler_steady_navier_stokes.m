function [H,b,Jac] = assembler_steady_navier_stokes(fespace_u,fespace_p,fun,nu,dirichlet_functions,neumann_functions)
% Assemble steady navier-stokes matrix and rhs with boundary conditions.
% The returned matrix is a function depending on the solution.
% input=
%           fespace_u: finite elemnet space for velocity
%           fespace_p: finite elemnet space for pressure
%           fun: anonymous function or scalar of the forcing term
%               If scalar, fun must be [0;0] and the code is optimized
%           nu: anonymous function or scalar for the diffusion term. If
%           scalar, the code is optimized on structured meshes
%           dirichlet_functions: Dirichlet boundary data
%           neumann_functions: Neumann_boundary data
% output=
%           H: system matrix
%           b: right handside
%           Jac: Jacobian of the matrix

[H,b] = assembler_steady_stokes(fespace_u,fespace_p,fun,nu,dirichlet_functions, ...
    neumann_functions);

H = @(u) add_convective_term(H,u,fespace_u);