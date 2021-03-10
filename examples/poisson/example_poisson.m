clear
clc
figure

% With this script, we solve the equation -mu * d_i^2u/d_{x_i}^2 = f.

bottom_left_corner_x = 0;
bottom_left_corner_y = 0;

% size of the domain (unitary square)
L = 1;
H = 1;

% number of elements in the two directions x and y
n_elements_x = 60;
n_elements_y = 60;

% create square mesh. To visualize it, type:
% draw_mesh(mesh)
mesh = create_mesh(bottom_left_corner_x, ...
                   bottom_left_corner_y, ...
                   L,H,n_elements_x,n_elements_y);

% boundary conditions. Each flag correspond to an edge of the square. We 
% start from the bottom and proceed counterclockwise. 
% 1: Dirichlet boundary condition (see dirichlet_functions) we impose the 
% value of the solution u
% 0: We prescribe Neumann conditions (see neumann_functions) we impose the 
% value of the normal derivative u.
% In this case, we are imposing homogeoneus ( = 0) Dirichlet boundary 
% conditions on the whole boundary.
bc_flags = [1 1 1 1];

% create fespace (i.e., the space of polynomial functions). In this case, 
% we create a space in which the local polynomial functions are cubic (P3).
% Other choices are: 
% 'P1': linear polynomial functions
% 'P2': quadratic polynomial functions
fespace = create_fespace(mesh,'P3',bc_flags);

% right hand side
f = @(x) sin(x(1,:).*x(2,:)).*x(2,:).^3;
mu = 3.5;
% 
% dirichlet_functions = @(x) [3*x(1);0;0;sin(x(2))];
% neumann_functions = @(x) [0;1;1;0];

dirichlet_functions = @(x) 0*[x(1);x(1);x(1);x(1)];
neumann_functions = @(x) [0;0;0;0];

% Here we assemble the elements of the system A*sol = b.
% A discretizes the Laplacian operator, b the right handside.
% sol is the the vector of degrees of freedom.
[A,b] = assembler_poisson(fespace,f,mu,dirichlet_functions,neumann_functions);

sol = A\b;

plot_fe_function(sol,fespace)
% export_vtk_scalar(sol,fespace,'example_poisson.vtk'); <- this is to
% export the function to Paraview. Not necessary for now.