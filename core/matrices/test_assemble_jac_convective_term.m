% common variables

tol = 1e-12;

%% Test 1: check if size is correct

mesh = create_mesh(0,0,1,1,20,20);
fespace_u = create_fespace(mesh,'P2',[0 0 0 0]);

fespace_u.mesh.type = '';

C1 = assemble_convective_term(fespace_u,ones(size(fespace_u.nodes,1)*2,1));
assert(size(C1,1) == size(fespace_u.nodes,1));

%% Test 2: check if optimized code is correct

mesh = create_mesh(0,0,1,1,20,21);
fespace_u = create_fespace(mesh,'P2',[0 0 0 0]);

vector_field = rand(size(fespace_u.nodes,1)*2,1);

C1 = assemble_jac_convective_term(fespace_u,vector_field,2,1);

fespace_u.mesh.type = '';

C2 = assemble_jac_convective_term(fespace_u,vector_field,2,1);

x = rand(size(C1,1),1);
assert(norm(C1*x-C2*x) < tol);

%% Test 3: check if Jacobian is equal to finite differences
epsilon = 1e-9;

mesh = create_mesh(0,0,1,1,6,7);
mesh.type = 'unstructured';
fespace_u = create_fespace(mesh,'P2',[0 0 0 0]);

u = rand(size(fespace_u.nodes,1)*2,1);

C = assemble_convective_term(fespace_u,u);

n = length(u);
zz = zeros(n,1);

matblock = [C zeros(n/2); zeros(n/2) C];

J_appr = zeros(n);
for i = 1:n
    zz(i) = epsilon;
    uz = u + zz;
    C_star = assemble_convective_term(fespace_u,uz);
    matblockstar = [C_star zeros(n/2); zeros(n/2) C_star];
    J_appr(:,i) = (matblockstar*uz-matblock*u)/epsilon;
    zz(i) = 0;
end

J11 = full(assemble_jac_convective_term(fespace_u,u,1,1));
J12 = full(assemble_jac_convective_term(fespace_u,u,1,2));
J21 = full(assemble_jac_convective_term(fespace_u,u,2,1));
J22 = full(assemble_jac_convective_term(fespace_u,u,2,2));

assert(norm(J11-J_appr(1:n/2,1:n/2)+C) < 1e-7)
assert(norm(J12-J_appr(1:n/2,n/2+1:end)) < 1e-7)
assert(norm(J21-J_appr(n/2+1:end,1:n/2)) < 1e-7)
assert(norm(J22-J_appr(n/2+1:end,n/2+1:end)+C) < 1e-7)

%% Test 3: check if Jacobian is equal to finite differences with Dir bc
epsilon = 1e-9;

mesh = create_mesh(0,0,1,1,3,3);
fespace_u = create_fespace(mesh,'P2',[1 1 1 1]);

u = rand(size(fespace_u.nodes,1)*2,1);

C = assemble_convective_term(fespace_u,u);

C = apply_dirichlet_bc_matrix(C,fespace_u,1);

n = length(u);
zz = zeros(n,1);

matblock = [C zeros(n/2); zeros(n/2) C];

J_appr = zeros(n);
for i = 1:n
    zz(i) = epsilon;
    uz = u + zz;
    C_star = assemble_convective_term(fespace_u,uz);
    C_star = apply_dirichlet_bc_matrix(C_star,fespace_u,1);

    matblockstar = [C_star zeros(n/2); zeros(n/2) C_star];
    J_appr(:,i) = (matblockstar*uz-matblock*u)/epsilon;
    zz(i) = 0;
end

J11 = full(apply_dirichlet_bc_matrix(assemble_jac_convective_term(fespace_u,u,1,1),fespace_u,0));
J12 = full(apply_dirichlet_bc_matrix(assemble_jac_convective_term(fespace_u,u,1,2),fespace_u,0));
J21 = full(apply_dirichlet_bc_matrix(assemble_jac_convective_term(fespace_u,u,2,1),fespace_u,0));
J22 = full(apply_dirichlet_bc_matrix(assemble_jac_convective_term(fespace_u,u,2,2),fespace_u,0));

assert(norm(J11-J_appr(1:n/2,1:n/2)+C) < 1e-7)
assert(norm(J12-J_appr(1:n/2,n/2+1:end)) < 1e-7)
assert(norm(J21-J_appr(n/2+1:end,1:n/2)) < 1e-7)
assert(norm(J22-J_appr(n/2+1:end,n/2+1:end)+C) < 1e-7)
