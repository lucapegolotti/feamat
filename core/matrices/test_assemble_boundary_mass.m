% common variables

tol = 1e-4;

%% Test 1: check if size is correct

mesh = create_mesh(0,0,1,1,20,20);
fespace = create_fespace(mesh,'P2',[0 0 0 0]);

A1 = assemble_boundary_mass(fespace,1);

assert(size(A1,1) == size(fespace.nodes,1))

%% Test 2: check if result is correct

f = @(x) sin(x(1));
mesh = read_mesh('square.msh');
fespace = create_fespace(mesh,'P2',[0 0 0 0]);
fun = project_function(fespace,f);

% check boundary 1
M = assemble_boundary_mass(fespace,1);
res = fun'*M*fun;
x = -0.5:0.01:0.5;
fx = sin(x).*sin(x);
realres = trapz(x,fx);
assert(abs(res - realres) < tol);

% check boundary 2
M = assemble_boundary_mass(fespace,2);
res = fun'*M*fun;
realres = sin(0.5)*sin(0.5);
assert(abs(res - realres) < tol);

