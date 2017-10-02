clear all
close all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 20;
n1 = 20;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);

f = @(t,x) [0;0];
center_nu = [L/4;H/2]';
nu = @(x) 1;
% dirichlet_functions = @(t,x) [0 0;0 0;0 0;0 sin(4*pi*t)]';
% neumann_functions = @(t,x) [0;0;0;0];
dirichlet_functions = @(t,x) [0 0;0 0;0 0;0 sin(8*t*pi)]';
% neumann_functions = @(t,x) [0 0;0 0;0 0;-cos(16*pi*t) 0]';
neumann_functions = @(t,x) [0 0;0 0;0 0;0 0]';


% Create finite element space
bc = [1 0 1 1]; 

poly_degree = 'P2';
fespace_u = create_fespace(mesh,poly_degree,bc);

poly_degree = 'P1';
fespace_p = create_fespace(mesh,poly_degree,bc);

sol = solver_unsteady_stokes(fespace_u,fespace_p,0,1,0.01,f,@(x) [0;0],nu,dirichlet_functions,neumann_functions);
%%
visualize_stokes_solution(sol,0.01)