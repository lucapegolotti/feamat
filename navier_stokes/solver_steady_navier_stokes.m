function [solution] = solver_steady_navier_stokes(fespace_u,fespace_p,fun,nu,dirichlet_functions,neumann_functions,varargin)
integrate_f = 1;
integrate_neumann = 1;
if (nargin >= 9)
    opts = varargin{1};
    integrate_f = opts.integrate_f;
    integrate_neumann = opts.integrate_neumann;
end

fun1 = @(x) fun(x)'*[1;0];
fun2 = @(x) fun(x)'*[0;1];
dir1 = @(x) dirichlet_functions(x)'*[1;0];
dir2 = @(x) dirichlet_functions(x)'*[0;1];
neu1 = @(x) neumann_functions(x)'*[1;0];
neu2 = @(x) neumann_functions(x)'*[0;1];

bc_flags_u = fespace_u.bc;
thereisneumann = 1;

n_nodes_u = size(fespace_u.nodes,1);
n_nodes_p = size(fespace_p.nodes,1);

nodes_u1 = 1:n_nodes_u;
nodes_u2 = n_nodes_u+1:n_nodes_u*2;

if (length(find(bc_flags_u)) == 4)
    thereisneumann = 0;
end

u1 = zeros(n_nodes_u,1);
u2 = zeros(n_nodes_u,1);
p = zeros(n_nodes_p,1);

u = [u1;u2;p];

% assemble constant matrices
A = assemble_stiffness(nu,fespace_u);
B1 = assemble_divergence(fespace_u,fespace_p,'dx');
B2 = assemble_divergence(fespace_u,fespace_p,'dy');

zero_mat_u = zeros(n_nodes_u);
zero_mat_p = zeros(n_nodes_p);

% A = apply_dirichlet_bc_matrix(A,fespace_u,1);
B1_u = B1';
B2_u = B2';
B1_u = apply_dirichlet_bc_matrix(B1_u,fespace_u,0);
B2_u = apply_dirichlet_bc_matrix(B2_u,fespace_u,0);

n_vertices = size(fespace_u.mesh.vertices,1);

b1 = zeros(n_nodes_u,1);
b2 = zeros(n_nodes_u,1);

if (integrate_f)
    b1 = assemble_rhs(fespace_u,fun1);
    b2 = assemble_rhs(fespace_u,fun2);
end

if (thereisneumann && integrate_neumann)
   b1 = apply_neumann_bc(fespace_u,b1,neu1); 
   b2 = apply_neumann_bc(fespace_u,b2,neu2); 
end

diagblock = @(u) apply_dirichlet_bc_matrix(A+assemble_convective_term(fespace_u,[u(nodes_u1);u(nodes_u2)]),fespace_u,1);

b1 = apply_dirichlet_bc_rhs(b1,fespace_u,dir1);
b2 = apply_dirichlet_bc_rhs(b2,fespace_u,dir2);

b = [b1;b2;zeros(n_nodes_p,1)];

% perform newtwon iterations until convergence
tol = 1e-8;
maxit = 1000;

evalres = @(u) [diagblock(u) zero_mat_u -B1_u; zero_mat_u diagblock(u) -B2_u; -B1 -B2 zero_mat_p]*u-b;

res = evalres(u);
err = norm(res);

count = 0;
while (err > tol && count < maxit)
    count = count+1;
    disp('----------------------------');
    disp(['Newton iteration ',num2str(count), ', residual norm = ',num2str(err)]);
    C = assemble_convective_term(fespace_u,[u1;u2]);
    J11 = assemble_Jac_convective_term(fespace_u,[u1;u2],1,1); 
    J12 = assemble_Jac_convective_term(fespace_u,[u1;u2],1,2); 
    J21 = assemble_Jac_convective_term(fespace_u,[u1;u2],2,1); 
    J22 = assemble_Jac_convective_term(fespace_u,[u1;u2],2,2);
    b11 = apply_dirichlet_bc_matrix(A+C+J11,fespace_u,1);
    b12 = apply_dirichlet_bc_matrix(J12,fespace_u,0);
    b21 = apply_dirichlet_bc_matrix(J21,fespace_u,0);
    b22 = apply_dirichlet_bc_matrix(A+C+J22,fespace_u,1);

    J = [b11 b12 -B1_u; b21 b22 -B2_u; -B1 -B2 zero_mat_p];
    u = u - J\res;
    res = evalres(u);
    err = norm(res);
end
disp(['Newton reached converged with norm of residual = ', num2str(err)]);

u1 = u(1:n_vertices);
u2 = u(n_nodes_u+1:n_nodes_u+1+n_vertices);

minnorm = min(sqrt(u1(1:n_vertices).^2+u2(1:n_vertices).^2));
maxnorm = max(sqrt(u1(1:n_vertices).^2+u2(1:n_vertices).^2));
minp = min(u(2*n_nodes_u+1:end));
maxp = max(u(2*n_nodes_u+1:end));

solution.u = u;
solution.u1 = u(1:n_nodes_u);
solution.u2 = u(n_nodes_u+1:n_nodes_u*2);
solution.p = u(n_nodes_u*2+1:end);
solution.mesh = fespace_u.mesh;
solution.fespace_u = fespace_u;
solution.fespace_p = fespace_p;
solution.minnorm = minnorm;
solution.maxnorm = maxnorm;
solution.minp = minp;
solution.maxp = maxp;
