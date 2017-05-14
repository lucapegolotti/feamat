function [H,b] = assembler_steady_stokes(fespace_u,fespace_p,fun,nu,dirichlet_functions,neumann_functions)
tic 
bc_flags_u = fespace_u.bc;
thereisneumann = 1;

n_nodes_u = size(fespace_u.nodes,1);
n_nodes_p = size(fespace_p.nodes,1);

if (length(find(bc_flags_u)) == 4)
    thereisneumann = 0;
end

fun1 = @(x) fun(x)'*[1;0];
fun2 = @(x) fun(x)'*[0;1];
dir1 = @(x) dirichlet_functions(x)'*[1;0];
dir2 = @(x) dirichlet_functions(x)'*[0;1];

A = assemble_stiffness(nu,fespace_u);
B1 = assemble_divergence(fespace_u,fespace_p,'dx');
B2 = assemble_divergence(fespace_u,fespace_p,'dy');
b1 = assemble_rhs(fespace_u,fun1);
b2 = assemble_rhs(fespace_u,fun2)*0;

zero_mat_u = zeros(n_nodes_u);
zero_mat_p = zeros(n_nodes_p);
zero_mat_up = zeros(n_nodes_u,n_nodes_p);

H1 = [A zero_mat_u -B1'];
H2 = [zero_mat_u A -B2'];
H3 = [-B1 -B2 zero_mat_p];
% H1 = [A zero_mat_u zero_mat_up];
% H2 = [zero_mat_u A zero_mat_up];
% H3 = [zero_mat_up' zero_mat_up' zero_mat_p];

[H1,b1] = apply_dirichlet_bc(H1,b1,fespace_u,dir1);
[H2(:,n_nodes_u+1:end),b2] = apply_dirichlet_bc(H2(:,n_nodes_u+1:end),b2,fespace_u,dir2);

if (thereisneumann)
   b1 = apply_neumann_bc(fespace_u,b1,neumann_functions); 
   b2 = apply_neumann_bc(fespace_u,b2,neumann_functions); 
end

H = [H1;H2;H3];
b = [b1;b2;zeros(n_nodes_p,1)];


elapsed = toc;
disp(['Assembly of the system took = ', num2str(elapsed),' s']);
disp('------------------------------');
