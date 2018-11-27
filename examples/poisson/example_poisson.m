clear all
clc
figure
bottom_left_corner_x = 0;
bottom_left_corner_y = 0;

L = 1;
H = 1;

n_elements_x = 60;
n_elements_y = 60;

mesh = create_mesh(bottom_left_corner_x, ...
                   bottom_left_corner_y, ...
                   L,H,n_elements_x,n_elements_y);

bc_flags = [1 1 1 1];
fespace = create_fespace(mesh,'P3',bc_flags);

f = @(x) sin(x(1,:).*x(2,:)).*x(2,:).^3;
mu = 3.5;
% 
% dirichlet_functions = @(x) [3*x(1);0;0;sin(x(2))];
% neumann_functions = @(x) [0;1;1;0];

% f = @(x) x(1,:)*0;
dirichlet_functions = @(x) 0*[x(1);x(1);x(1);x(1)];
neumann_functions = @(x) [0;0;0;0];

[A,b] = assembler_poisson(fespace,f,mu,dirichlet_functions,neumann_functions);

sol = A\b;

plot_fe_function(sol,fespace)
% export_vtk_scalar(sol,fespace,'example_poisson.vtk');