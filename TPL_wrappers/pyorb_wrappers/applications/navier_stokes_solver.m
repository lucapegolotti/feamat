function [sol] = navier_stokes_solver( param, fem_specifics )

mesh_file = fem_specifics.mesh_name;

mesh = read_mesh( mesh_file );

bc_flags    = [1 0 1 0 1 1];


fespace_u = create_fespace(mesh,'P2',bc_flags);
fespace_p = create_fespace(mesh,'P1',bc_flags);

f = [0;0];
mu = 0.1;

U = 1.0 * param(1);
r = 0.4;

v_in = @(x) (-x(2).^2+r^2)*U*4;

dirichlet_stokes_functions    = @(x) [0 0; 0 0; 0 0; 0 0; 0 0; v_in(x) 0]';
neumann_functions             = @(x) [0 0;0 0;0 0;0 0;0 0;0 0]';

[A_stokes, b_stokes] = assembler_steady_stokes( fespace_u, fespace_p, f, mu, dirichlet_stokes_functions, ...
                                         neumann_functions);
u_lifting_stokes = A_stokes \ b_stokes;



% n_nodes_u = size(fespace_u.nodes,1);
% u_lifting_stokes = u_lifting_stokes(1:2*n_nodes_u);

dirichlet_functions = @(x) [0 0; 0 0; 0 0; 0 0; 0 0; 0 0]';

C_1 = A_stokes * 0.0;
C_1 = add_convective_term( C_1, u_lifting_stokes, fespace_u );
C_1 = add_flipped_convective_term( C_1, u_lifting_stokes, fespace_u );

[A_no_lifting, b_no_lifting] = assembler_steady_navier_stokes( fespace_u, fespace_p, f, mu, dirichlet_functions,...
                                                               neumann_functions );
n_nodes_u = size(fespace_u.nodes,1);

% building rhs of the equation coming from lifting
b = b_no_lifting - A_no_lifting(u_lifting_stokes) * u_lifting_stokes;
b(1:n_nodes_u) = apply_dirichlet_bc_rhs( b(1:n_nodes_u), fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
b(n_nodes_u+1:2*n_nodes_u) = apply_dirichlet_bc_rhs( b(n_nodes_u+1:2*n_nodes_u), fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );

A = @(u) A_no_lifting(u) + C_1;


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

sol_lifting = sol;
sol_lifting.u1 = sol.u1 + u_lifting_stokes(1:n_nodes_u);
sol_lifting.u2 = sol.u2 + u_lifting_stokes(n_nodes_u+1:2*n_nodes_u);
plot_fe_fluid_function(sol_lifting,'U','contourf');
% export_vtk_fluid(sol,'test','U')
axis equal 

norm(sol_lifting.u1)
norm(sol_lifting.u2)


end