clear all
close all
clc

% Set dimension of the domain and parameters of the mesh
L = 4;
H = 1;

n2 = 20;
n1 = 80;

% Create and display the mesh
mesh = create_mesh(L,H,n1,n2);

f = @(t,x) [0;0];
center_nu = [L/4;H/2]';
nu = @(x) 1;
dirichlet_functions = @(t,x) [0 0;0 0;0 0;0 0]';
neumann_functions = @(t,x) [0 0;0 0;0 0;cos(8*pi*t)*10 0]';

% Create finite element space
bc = [1 0 1 0]; 

poly_degree = 'P2';
fespace_u = create_fespace(mesh,poly_degree,bc);

poly_degree = 'P1';
fespace_p = create_fespace(mesh,poly_degree,bc);

% options
opts.integrate_f = 0;
opts.integrate_neumann = 1;

sol = solver_navier_stokes(fespace_u,fespace_p,0,1,0.01,f,@(x) [0;0],nu,dirichlet_functions,neumann_functions,opts);
%%
visualize_stokes_solution(sol,0.01)