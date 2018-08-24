function [H,b] = assembler_steady_stokes(fespace_u,fespace_p,fun,nu,dirichlet_functions,neumann_functions)
% Assemble steady stokes matrix and rhs with boundary conditions
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


bc_flags_u = fespace_u.bc;
thereisneumann = 1;

n_nodes_u = size(fespace_u.nodes,1);
n_nodes_p = size(fespace_p.nodes,1);

if (length(find(bc_flags_u)) == length(bc_flags_u))
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

A = assemble_stiffness(nu,fespace_u);
B1 = assemble_divergence(fespace_u,fespace_p,'dx');
B2 = assemble_divergence(fespace_u,fespace_p,'dy');

if (integratef)
    b1 = assemble_rhs(fespace_u,fun1);
    b2 = assemble_rhs(fespace_u,fun2);
else
    b1 = zeros(n_nodes_u,1);
    b2 = zeros(n_nodes_u,1);
end

zero_mat_u = sparse(n_nodes_u,n_nodes_u);
zero_mat_p = sparse(n_nodes_p,n_nodes_p);
zero_mat_up = sparse(n_nodes_u,n_nodes_p);

H1 = [A zero_mat_u' -B1'];
H2 = [zero_mat_u A -B2'];
H3 = [-B1 -B2 zero_mat_p];

if (thereisneumann)
    b1 = apply_neumann_bc(b1,fespace_u,neu1);
    b2 = apply_neumann_bc(b2,fespace_u,neu2);
else
    % arbitrarly fix pressure in first dof if there is no Neumann condition
    H3(1,:) = 0;
    H3(1,2*n_nodes_u+1) = 1;
end

[H1,b1] = apply_dirichlet_bc(H1,b1,fespace_u,dir1);
[H2(:,n_nodes_u+1:end),b2] = apply_dirichlet_bc(H2(:,n_nodes_u+1:end),b2,fespace_u,dir2);

H = [H1;H2;H3];
b = [b1;b2;zeros(n_nodes_p,1)];