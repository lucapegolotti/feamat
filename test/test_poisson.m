clear all
close all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 20;
n1 = n2*L;

% Create and display the mesh
mesh = create_mesh(L,H,n1,n2);
draw_mesh(mesh);

f = @(x,y) 2*pi^2*sin(pi*x)*sin(pi*y);
mu = @(x,y) 1;
dirichlet_functions = @(x,y) [0;0;0;0];

% Create finite element space
bc = [1 1 1 1]; 

poly_degree = 'P1';
fespace = create_fespace(mesh,poly_degree,bc);

[A,b] = assembler_poisson(f,mu,fespace,dirichlet_functions);


figure()
spy(A)

sol = A\b;

n1 = size(mesh.X,1);
n2 = size(mesh.X,2);


surf(mesh.X,mesh.Y,reshape(sol,n1,n2));

% mat = [-3 3; 3 0];
% grad = mat*[0;1];
% hold on
% plot([0 grad(2)]+1/3,[0 grad(1)]+1/3)


