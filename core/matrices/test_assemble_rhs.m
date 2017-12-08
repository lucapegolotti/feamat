% common variables

tol = 1e-12;

%% Test 1: check if size is correct

mesh = create_mesh(0,0,1,1,20,20);
fespace = create_fespace(mesh,'P2',[0 0 0 0]);

rhs1 = assemble_rhs(fespace, 0);
rhs2 = assemble_rhs(fespace, @(x) x(1)*x(2));

assert(length(rhs1) == size(fespace.nodes,1))
assert(length(rhs2) == size(fespace.nodes,1))
