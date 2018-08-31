function [A,b,uLift,iNodes] = assembler_poisson_lifting(fespace,fun,mu,dirichlet_functions,neumann_functions)
% Assemble poisson matrix with lifting boundary conditions
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
%           uLift: lifting function
%           iNodes: internal Nodes
%
%   Author: Stefano Pagani <stefano.pagani at polimi.it>

bc_flags = fespace.bc;

thereisneumann = 1;

if (length(find(bc_flags)) == length(bc_flags))
    thereisneumann = 0;
end

A = assemble_stiffness(mu,fespace);
b = assemble_rhs(fespace,fun);

if (thereisneumann)
   b = apply_neumann_bc(b,fespace,neumann_functions); 
end

% SP: apply Dirichlet boundary conditions using lifting 
[A,b,uLift,iNodes] = apply_dirichlet_bc_lifting(A,b,fespace,dirichlet_functions);