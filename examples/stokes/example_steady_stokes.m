clear all
clc

bottom_left_corner_x = 0;
bottom_left_corner_y = 0;

L = 1;
H = 1;

n_elements_x = 100;
n_elements_y = 100;

mesh = create_mesh(bottom_left_corner_x, ...
                   bottom_left_corner_y, ...
                   L,H,n_elements_x,n_elements_y);

bc_flags = [1 1 1 1];

fespace_u = create_fespace(mesh,'P2',bc_flags);
fespace_p = create_fespace(mesh,'P1',bc_flags);

f = [0;0];
mu = 1;

dirichlet_functions = @(x) [0 0;0 0;1 0;0 0]';
neumann_functions = @(x) [0 0;0 0;0 0;0 0]';

[A,b] = assembler_steady_stokes(fespace_u,fespace_p,f,mu,dirichlet_functions,neumann_functions);

sol = solve_fluid_system(A,b,fespace_u,fespace_p);

% we plot the result using matlab and by exporting it to vtk
plot_fe_fluid_function(sol,'U');
export_vtk_fluid(sol,'example_steady_stokes.vtk')