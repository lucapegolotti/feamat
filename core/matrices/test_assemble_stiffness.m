% common variables

tol = 1e-12;

%% Test 1: check if size is correct

mesh = create_mesh(0,0,1,1,20,20);
fespace = create_fespace(mesh,'P2',[0 0 0 0]);

A1 = assemble_stiffness(1, fespace);

assert(size(A1,1) == size(fespace.nodes,1))

%% Test 2: check if optimized code is correct

mesh = create_mesh(0,0,1,1,20,25);
fespace = create_fespace(mesh,'P2',[0 0 0 0]);

A1 = assemble_stiffness(1, fespace);

A2 = assemble_stiffness(@(x) 1,fespace);

fespace.mesh.type = '';
A3 = assemble_stiffness(@(x) 1, fespace);

x = rand(size(A1,1),1);

assert(norm(A1*x-A2*x) < tol);
assert(norm(A1*x-A3*x) < tol);

