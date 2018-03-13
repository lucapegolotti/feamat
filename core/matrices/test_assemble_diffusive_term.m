% test diffusive term assembly
% common variables
tol = 1e-12;
%% first test: dimensional consistency
mesh = create_mesh(0,0,1,1,20,20);
fespace = create_fespace(mesh,'P2',[0 0 0 0]);
mu = eye(2);
A1 = assemble_diffusive_term(fespace, mu);
assert(size(A1,1) == size(fespace.nodes,1))
%% second test: isotropic unitary diffusivity
mesh = create_mesh(0,0,1,1,20,20);
fespace = create_fespace(mesh,'P2',[0 0 0 0]);
mu = eye(2);
A1 = assemble_diffusive_term(fespace, mu);
A2 = assemble_stiffness(1, fespace);
assert(norm(A1-A2,Inf) < tol)
%% third test: isotropic diffusion in rectangle
mesh = create_mesh(0,0,0.5,1,10,20);
fespace = create_fespace(mesh,'P2',[0 0 0 0]);
mu = eye(2);
A1 = assemble_diffusive_term(fespace, mu);
A2 = assemble_stiffness(1, fespace);
assert(norm(A1-A2,Inf) < tol)
%% fourth test: check constant and anonymous function yield same result
mesh = create_mesh(0,0,1,1,10,10);
fespace = create_fespace(mesh,'P2',[0 0 0 0]);
mu_const = eye(2);
mu_fun = @(x) [1 0;0 1];
A1 = assemble_diffusive_term(fespace, mu_const);
A2 = assemble_diffusive_term(fespace, mu_fun);
% use machine epsilon as tolerance, since the code is identical
assert(norm(A1-A2,Inf) < eps)
%% fifth test: check convergence order isotropic case
N = [5 10];
errH1 = zeros(size(N));
flags = [1 1 1 1];
dir_functions = @(x) [0;0;0;0];
mu = eye(2);
uexxy = @(x,y) sin(pi*x).*sin(pi*y);
duexdx = @(x,y) pi*cos(pi*x).*sin(pi*y);
duexdy = @(x,y) pi*sin(pi*x).*cos(pi*y);
uex = @(x) uexxy(x(1,:),x(2,:));
graduex = @(x) [duexdx(x(1,:),x(2,:));duexdy(x(1,:),x(2,:))];
funxy = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
fun = @(x) funxy(x(1,:),x(2,:));
for i = 1:length(N)
    mesh = create_mesh(0,0,1,1,N(i),N(i));
    fespace = create_fespace(mesh,'P2',flags);
    A = assemble_diffusive_term(fespace, mu);
    f = assemble_rhs(fespace, fun);
    A = apply_dirichlet_bc_matrix(A, fespace, 1);
    f = apply_dirichlet_bc_rhs(f, fespace, dir_functions);
    u = A \ f;
    errH1(i) = compute_H1_error(fespace, u, uex, graduex);
end
convH1 = log2(errH1(1)./errH1(2));
% check relative error on the convergence order
assert(abs(convH1 - 2.0) / 2.0 < 0.01);
%% sixth test: check convergence order anisotropic uniform case
N = [5 10];
errH1 = zeros(size(N));
flags = [1 1 1 1];
dir_functions = @(x) [0;0;0;0];
alpha = -0.3;
mu = [1 alpha; alpha 1];
uexxy = @(x,y) sin(pi*x).*sin(pi*y);
duexdx = @(x,y) pi*cos(pi*x).*sin(pi*y);
duexdy = @(x,y) pi*sin(pi*x).*cos(pi*y);
uex = @(x) uexxy(x(1,:),x(2,:));
graduex = @(x) [duexdx(x(1,:),x(2,:));duexdy(x(1,:),x(2,:))];
funxy = @(x,y) pi^2 * (2*sin(pi*x).*sin(pi*y) - 2 * alpha * cos(pi*x).*cos(pi*y));
fun = @(x) funxy(x(1,:),x(2,:));
for i = 1:length(N)
    mesh = create_mesh(0,0,1,1,N(i),N(i));
    fespace = create_fespace(mesh,'P2',flags);
    A = assemble_diffusive_term(fespace, mu);
    f = assemble_rhs(fespace, fun);
    A = apply_dirichlet_bc_matrix(A, fespace, 1);
    f = apply_dirichlet_bc_rhs(f, fespace, dir_functions);
    u = A \ f;
    errH1(i) = compute_H1_error(fespace, u, uex, graduex);
end
convH1 = log2(errH1(1)./errH1(2));
% check relative error on the convergence order
assert(abs(convH1 - 2.0) / 2.0 < 0.01);
%% seventh test: check convergence order anisotropic non-uniform case
N = [5 10];
errH1 = zeros(size(N));
flags = [1 1 1 1];
dir_functions = @(x) [0;0;0;0];
alpha = 1.0;
mu11 = @(x,y) 1 + x;
mu12 = @(x,y) alpha * x .* y;
mu22 = @(x,y) 1 + y;
mu = @(x) [mu11(x(1,:),x(2,:)) mu12(x(1,:),x(2,:)); mu12(x(1,:),x(2,:)) mu22(x(1,:),x(2,:))];
uexxy = @(x,y) sin(pi*x) .* sin(pi*y);
duexdx = @(x,y) pi * cos(pi*x) .* sin(pi*y);
duexdy = @(x,y) pi * sin(pi*x) .* cos(pi*y);
uex = @(x) uexxy(x(1,:),x(2,:));
graduex = @(x) [duexdx(x(1,:),x(2,:));duexdy(x(1,:),x(2,:))];
funxy = @(x,y) pi^2 * (2 + x + y) .* sin(pi*x) .* sin(pi*y) + ...
               -2 * alpha * pi^2 * x .* y .* cos(pi*x) .* cos(pi*y) +...
               -(1 + alpha*x) .* pi .* cos(pi*x) .* sin(pi*y) + ...
               -(1 + alpha*y) .* pi .* sin(pi*x) .* cos(pi*y);
fun = @(x) funxy(x(1,:),x(2,:));
for i = 1:length(N)
    mesh = create_mesh(0,0,1,1,N(i),N(i));
    fespace = create_fespace(mesh,'P2',flags);
    A = assemble_diffusive_term(fespace, mu);
    f = assemble_rhs(fespace, fun);
    A = apply_dirichlet_bc_matrix(A, fespace, 1);
    f = apply_dirichlet_bc_rhs(f, fespace, dir_functions);
    u = A \ f;
    errH1(i) = compute_H1_error(fespace, u, uex, graduex);
end
convH1 = log2(errH1(1)/errH1(2));
% check relative error on the convergence order
assert(abs(convH1 - 2.0) / 2.0 < 0.01);
%% eighth test: check convergence order anisotropic uniform case non-structured mesh
N = [5 10];
errH1 = zeros(size(N));
flags = [1 1 1 1];
dir_functions = @(x) [0;0;0;0];
alpha = 0.6;
mu = [1 alpha; alpha 1];
uexxy = @(x,y) sin(pi*x).*sin(pi*y);
duexdx = @(x,y) pi*cos(pi*x).*sin(pi*y);
duexdy = @(x,y) pi*sin(pi*x).*cos(pi*y);
uex = @(x) uexxy(x(1,:),x(2,:));
graduex = @(x) [duexdx(x(1,:),x(2,:));duexdy(x(1,:),x(2,:))];
funxy = @(x,y) pi^2 * (2*sin(pi*x).*sin(pi*y) - 2 * alpha * cos(pi*x).*cos(pi*y));
fun = @(x) funxy(x(1,:),x(2,:));
for i = 1:length(N)
    mesh = create_mesh(0,0,1,1,N(i),N(i));
    mesh.type = '';
    fespace = create_fespace(mesh,'P2',flags);
    A = assemble_diffusive_term(fespace, mu);
    f = assemble_rhs(fespace, fun);
    A = apply_dirichlet_bc_matrix(A, fespace, 1);
    f = apply_dirichlet_bc_rhs(f, fespace, dir_functions);
    u = A \ f;
    errH1(i) = compute_H1_error(fespace, u, uex, graduex);
end
convH1 = log2(errH1(1)./errH1(2));
% check relative error on the convergence order
assert(abs(convH1 - 2.0) / 2.0 < 0.01);
%% ninth test: check convergence order anisotropic non-uniform case non-structured mesh
N = [5 10];
errH1 = zeros(size(N));
flags = [1 1 1 1];
dir_functions = @(x) [0;0;0;0];
alpha = 1.0;
mu11 = @(x,y) 1 + x;
mu12 = @(x,y) alpha * x .* y;
mu22 = @(x,y) 1 + y;
mu = @(x) [mu11(x(1,:),x(2,:)) mu12(x(1,:),x(2,:)); mu12(x(1,:),x(2,:)) mu22(x(1,:),x(2,:))];
uexxy = @(x,y) sin(pi*x) .* sin(pi*y);
duexdx = @(x,y) pi * cos(pi*x) .* sin(pi*y);
duexdy = @(x,y) pi * sin(pi*x) .* cos(pi*y);
uex = @(x) uexxy(x(1,:),x(2,:));
graduex = @(x) [duexdx(x(1,:),x(2,:));duexdy(x(1,:),x(2,:))];
funxy = @(x,y) pi^2 * (2 + x + y) .* sin(pi*x) .* sin(pi*y) + ...
               -2 * alpha * pi^2 * x .* y .* cos(pi*x) .* cos(pi*y) +...
               -(1 + alpha*x) .* pi .* cos(pi*x) .* sin(pi*y) + ...
               -(1 + alpha*y) .* pi .* sin(pi*x) .* cos(pi*y);
fun = @(x) funxy(x(1,:),x(2,:));
for i = 1:length(N)
    mesh = create_mesh(0,0,1,1,N(i),N(i));
    mesh.type = '';
    fespace = create_fespace(mesh,'P2',flags);
    A = assemble_diffusive_term(fespace, mu);
    f = assemble_rhs(fespace, fun);
    A = apply_dirichlet_bc_matrix(A, fespace, 1);
    f = apply_dirichlet_bc_rhs(f, fespace, dir_functions);
    u = A \ f;
    errH1(i) = compute_H1_error(fespace, u, uex, graduex);
end
convH1 = log2(errH1(1)./errH1(2));
% check relative error on the convergence order
assert(abs(convH1 - 2.0) / 2.0 < 0.01);