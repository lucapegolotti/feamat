clear all
close all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 30;
n1 = 30;

% Create and display the mesh
mesh = create_mesh(L,H,n1,n2);
% 

f = @(x) [0;0];
%nu = @(x) 1.*(norm(x-[1;0.5])>0.2) + 1e-14.*(norm(x-[1;0.5])<=0.2);
nu = @(x) 1;
dirichlet_functions = @(x) [0 0;0 0;0 0;0 1]';
neumann_functions = @(x) [0;1;0;0];

% Create finite element space
bc = [1 0 0 1]; 

poly_degree = 'P2';
fespace_u = create_fespace(mesh,poly_degree,bc);

poly_degree = 'P1';
fespace_p = create_fespace(mesh,poly_degree,bc);

% Assemble matrix and rhs
[A,b] = assembler_steady_stokes(fespace_u, fespace_p, f, nu, dirichlet_functions, neumann_functions);

% Solve the linear system
% A(1,:) = zeros(1,size(A,2));
% A(1,1) = 1;
% b(1) = 1;

sol = A\b;


n_nodes_u = size(fespace_u.nodes,1);

u1 = sol(1:n_nodes_u);
u2 = sol(n_nodes_u+1:2*n_nodes_u);
p = sol(n_nodes_u*2+1:end);

subplot(1,2,1)
plot_solution_vp(fespace_u,fespace_p,sol,'U')
title('V')
axis([0 L 0 H])
%hold on
%draw_mesh(mesh);

pbaspect([L H 1])
subplot(1,2,2)
plot_solution_vp(fespace_u,fespace_p,sol,'P')
%hold on
%draw_mesh(mesh);
title('P')
pbaspect([L H 1])

% l2norm = compute_norm(fespace,sol,'L2');
% display(['Norm = ', num2str(l2norm)]);