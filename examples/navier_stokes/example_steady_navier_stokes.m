clear all
clc

bottom_left_corner_x = 0;
bottom_left_corner_y = 0;

L = 1;
H = 1;

n_elements_x = 40;
n_elements_y = 40;

mesh = create_mesh(bottom_left_corner_x, ...
    bottom_left_corner_y, ...
    L,H,n_elements_x,n_elements_y);

bc_flags = [1 1 1 1];

fespace_u = create_fespace(mesh,'P2',bc_flags);
fespace_p = create_fespace(mesh,'P1',bc_flags);

f = [0;0];
mu = 1;

U = 200;

dirichlet_functions = @(x) [0 0;0 0;U 0;0 0]';
neumann_functions = @(x) [0 0;0 0;0 0;0 0]';

[A,b] = assembler_steady_navier_stokes(fespace_u,fespace_p,f,mu,dirichlet_functions,...
    neumann_functions);

% solve stokes problem with u = 0 to get initial guess for newton's method
u0 = zeros(size(fespace_u.nodes,1)*2,1);
x0 = A(u0)\b;

% solve system with newton's method
method.name = 'newton';
method.f = @(u) A(u)*u-b;
method.x0 = x0;
method.jac = @(u) build_jac_navier_stokes(A,u,fespace_u);
method.tol = 1e-8;
method.maxit = 100;

[sol,err,it] = solve_fluid_system(A,b,fespace_u,fespace_p,method);

% plot and export
plot_fe_fluid_function(sol,'U');
%export_vtk_fluid(sol,'example_steady_navier_stokes.vtk');