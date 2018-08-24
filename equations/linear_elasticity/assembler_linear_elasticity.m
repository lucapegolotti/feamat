function [H,b] = assembler_linear_elasticity(fespace,fun,poisson,young,dirichlet_functions,neumann_functions)
% Assemble steady stokes matrix and rhs with boundary conditions
% input=
%           fespace: finite elemnet space for displacement
%           fun: anonymous function or scalar of the forcing term
%           poisson: poisson ration
%           young: young modulus
%           dirichlet_functions: Dirichlet boundary data
%           neumann_functions: Neumann_boundary data
% output=
%           H: system matrix
%           b: right handside


if (strcmp(fespace.mesh.type,'structured'))
    warning(['Attention! Assemble of structure matrices is not optimized for', ...
             ' structured meshes']);
end

lambda = poisson * young/((1 - 2*poisson)*(1+poisson));
mu = young/(2*(1+poisson));

bc_flags = fespace.bc;
thereisneumann = 1;

n_nodes = size(fespace.nodes,1);

if (length(find(bc_flags)) == length(bc_flags))
    thereisneumann = 0;
end

integratef = 1;
if (~isa(fun,'function_handle'))
    if (fun(1) == 0 && fun(1) == 0)
        integratef = 0;
    else
        error('Constant fun can only be zero! Pass it as anonymous function instead');
    end
end

fun1 = @(x) [1 0]*fun(x);
fun2 = @(x) [0 1]*fun(x);
dir1 = @(x) dirichlet_functions(x)'*[1;0];
dir2 = @(x) dirichlet_functions(x)'*[0;1];
neu1 = @(x) neumann_functions(x)'*[1;0];
neu2 = @(x) neumann_functions(x)'*[0;1];

A11 = assemble_structure_single_matrix(fespace,1,1);
A12 = assemble_structure_single_matrix(fespace,1,2);
A21 = assemble_structure_single_matrix(fespace,2,1);
A22 = assemble_structure_single_matrix(fespace,2,2);

if (integratef)
    b1 = assemble_rhs(fespace,fun1);
    b2 = assemble_rhs(fespace,fun2);
else
    b1 = zeros(n_nodes,1);
    b2 = zeros(n_nodes,1);
end

A_11 = (lambda + 2*mu)*A11 + mu * A22;
A_11= apply_dirichlet_bc_matrix(A_11,fespace,1);

A_12 = lambda*A21 + mu*A12;
A_12  = apply_dirichlet_bc_matrix(A_12,fespace,0);

A_21 = lambda*A12 + mu*A21;
A_21  = apply_dirichlet_bc_matrix(A_21,fespace,0);

A_22 = (lambda + 2*mu)*A22 + mu * A11;
A_22= apply_dirichlet_bc_matrix(A_22,fespace,1);

H = [A_11 A_12;
     A_21 A_22];
  
if (thereisneumann)
    b1 = apply_neumann_bc(b1,fespace,neu1);
    b2 = apply_neumann_bc(b2,fespace,neu2);
end

b1 = apply_dirichlet_bc_rhs(b1,fespace,dir1);
b2 = apply_dirichlet_bc_rhs(b2,fespace,dir2);
  
b = [b1;b2];