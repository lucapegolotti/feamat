% common variables
tol = 1e-12;

%% Test 1: check rhs stays equal with homogeneous data

mesh = create_mesh(0,0,1,1,10,10);
fespace = create_fespace(mesh,'P2',[0 1 0 0]);
n_nodes = size(fespace.nodes,1);

b = ones(n_nodes,1);
b = apply_neumann_bc(b,fespace,@(x) [0;0;0;0]);

assert(norm(b-b) < tol)

%% Test2: check rhs is integrated correctly with P1 and constant function

mesh = create_mesh(0,0,1,1,1,1);
fespace = create_fespace(mesh,'P1',[1 0 1 1]);
n_nodes = size(fespace.nodes,1);

b = zeros(n_nodes,1);
b = apply_neumann_bc(b,fespace,@(x) [0 1 0 0]);

assert(b(2) == 0.5);
assert(b(4) == 0.5);

%% Test2: check rhs is integrated correctly with P2 and constant function

mesh = create_mesh(0,0,1,1,1,1);
fespace = create_fespace(mesh,'P2',[1 0 1 1]);
n_nodes = size(fespace.nodes,1);

b = zeros(n_nodes,1);
b = apply_neumann_bc(b,fespace,@(x) [0 1 0 0]);

assert(abs(b(2) - 1/6) < tol);
assert(abs(b(4) - 1/6) < tol);
assert(abs(b(6)/4 - 1/6) < tol);


