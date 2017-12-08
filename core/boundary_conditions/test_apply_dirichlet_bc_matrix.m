% common variables

%% Test 1: check if matrix is diagonalized with 2 elements

mesh = create_mesh(0,0,1,1,1,1);
fespace = create_fespace(mesh,'P1',[1 0 0 0]);
n_nodes = size(fespace.nodes,1);

A = rand(n_nodes);
A = apply_dirichlet_bc_matrix(A,fespace,1);

diag_index = 1;
assert(A(diag_index,diag_index) == 1);
assert(sum(A(diag_index,:)) == 1)
diag_index = 2;
assert(A(diag_index,diag_index) == 1);
assert(sum(A(diag_index,:)) == 1)

fespace = create_fespace(mesh,'P1',[0 1 0 0]);

A = rand(n_nodes);
A = apply_dirichlet_bc_matrix(A,fespace,1);

diag_index = 2;
assert(A(diag_index,diag_index) == 1);
assert(sum(A(diag_index,:)) == 1)
diag_index = 4;
assert(A(diag_index,diag_index) == 1);
assert(sum(A(diag_index,:)) == 1)

fespace = create_fespace(mesh,'P1',[0 0 1 0]);

A = rand(n_nodes);
A = apply_dirichlet_bc_matrix(A,fespace,1);

diag_index = 3;
assert(A(diag_index,diag_index) == 1);
assert(sum(A(diag_index,:)) == 1)
diag_index = 4;
assert(A(diag_index,diag_index) == 1);
assert(sum(A(diag_index,:)) == 1)

fespace = create_fespace(mesh,'P1',[0 0 0 1]);

A = rand(n_nodes);
A = apply_dirichlet_bc_matrix(A,fespace,1);

diag_index = 1;
assert(A(diag_index,diag_index) == 1);
assert(sum(A(diag_index,:)) == 1)
diag_index = 3;
assert(A(diag_index,diag_index) == 1);
assert(sum(A(diag_index,:)) == 1)
