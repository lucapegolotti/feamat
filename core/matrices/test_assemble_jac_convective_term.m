% common variables

tol = 1e-12;

%% Test 1: check if size is correct

mesh = create_mesh(0,0,1,1,20,20);
fespace_u = create_fespace(mesh,'P2',[0 0 0 0]);

fespace_u.mesh.type = '';

C1 = assemble_convective_term(fespace_u,ones(size(fespace_u.nodes,1)*2,1));
assert(size(C1,1) == size(fespace_u.nodes,1));

%% Test 2: check if optimized code is correct

mesh = create_mesh(0,0,1,1,20,20);
fespace_u = create_fespace(mesh,'P2',[0 0 0 0]);

vector_field = rand(size(fespace_u.nodes,1)*2,1);

C1 = assemble_jac_convective_term(fespace_u,vector_field,2,1);

fespace_u.mesh.type = '';

C2 = assemble_jac_convective_term(fespace_u,vector_field,2,1);

x = rand(size(C1,1),1);
assert(norm(C1*x-C2*x) < tol);
