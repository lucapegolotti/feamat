% common variables

tol = 1e-12;

%% Test 1: check if size is correct

mesh = create_mesh(0,0,1,1,20,20);
fespace_u = create_fespace(mesh,'P2',[0 0 0 0]);
fespace_p = create_fespace(mesh,'P1',[0 0 0 0]);

B = assemble_divergence(fespace_u,fespace_p,'dx');

assert(size(B,1) == size(fespace_p.nodes,1))
assert(size(B,2) == size(fespace_u.nodes,1))

%% Test 2: check if optimized code is correct

mesh = create_mesh(0,0,1,1,20,20);
fespace_u = create_fespace(mesh,'P2',[0 0 0 0]);
fespace_p = create_fespace(mesh,'P1',[0 0 0 0]);

B1 = assemble_divergence(fespace_u,fespace_p,'dx');

fespace_u.mesh.type  = '';
B2 = assemble_divergence(fespace_u,fespace_p,'dx');

x = rand(size(fespace_u.nodes,1),1);
assert(norm(B1*x - B2*x) < tol)
