% clear all
% close all
% clc

% Set dimension of the domain and parameters of the mesh
L = 3;
H = 3;

n2 = 50;
n1 = n2;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);
draw_mesh(mesh);

f = @(x) sin(x(1))*cos(x(2));
mu = @(x) 1;
dirichlet_functions = @(x) [0;0;0;x(2)*(1-x(2))];
neumann_functions = @(x) [0;0;0;0];

% Create finite element space
bc = [0 0 0 0]; 

poly_degree = 'P1';
fespace = create_fespace(mesh,poly_degree,bc);

% Assemble matrix and rhs
[A,b] =   assembler_poisson(fespace,f,mu,dirichlet_functions,neumann_functions);

v = assemble_vector_lagrange_multiplier(fespace);
mat = [A v; v' 0];

% Solve the linear system
sol = mat\[b;10];

lmult = sol(end);

sol = sol(1:end-1);

figure
plot_solution_on_fespace(fespace,sol)
pbaspect([1 1 1])

integral = compute_integral_over_mesh(sol,fespace);
display(['Integral = ', num2str(integral)]);

integral = compute_integral_over_mesh(sol,fespace);
display(['Lagrange multiplier = ', num2str(lmult)]);