% common variables

tol = 1e-2;

%% Test 1: all fields assigned correctly, without physical entities
boundary_indicators = @(x) [(norm(x) < 0.3 + tol);
                            abs(x(1)-0.5) < tol;
                            abs(x(1)+0.5) < tol;
                            abs(x(2)-0.5) < tol;
                            abs(x(2)+0.5) < tol];

mesh = read_mesh('square_hole.msh',boundary_indicators);

assert(mesh.xp == -0.5);
assert(mesh.yp == -0.5); 
assert(mesh.L == 1);
assert(mesh.H == 1);
assert(abs(mesh.h-0.1086) < tol)

%% Test 2: all fields assigned correctly, with physical entities
boundary_indicators = @(x) [(norm(x) < 0.3 + tol);
                            abs(x(1)-0.5) < tol;
                            abs(x(1)+0.5) < tol;
                            abs(x(2)-0.5) < tol;
                            abs(x(2)+0.5) < tol];

mesh = read_mesh('bifurcation.msh',boundary_indicators);

assert(mesh.xp == 0);
assert(mesh.yp == -0.9); 
assert(mesh.L == 2.2);
assert(mesh.H == 1.8);
assert(abs(mesh.h-0.0589) < tol)

%% Test 3: matrices with structured and unstructured meshes are similar

xp = -0.5;
yp = -0.5;
L = 1;
H = 1;

test_function = @(x) sin(x(1)) .* cos(x(2));

mesh_str = create_mesh(xp,yp,L,H,30,30);
mesh_unstr = read_mesh('square.msh');

% structured fespace
fespace_str = create_fespace(mesh_str,'P2',[1 1 1 1]);
test_str = project_function(fespace_str,test_function);

fespace_unstr = create_fespace(mesh_unstr,'P2',[1 1 1 1]);
test_unstr = project_function(fespace_unstr,test_function);

A_str = assemble_stiffness(1,fespace_str);
A_unstr = assemble_stiffness(1,fespace_unstr);

res_str = A_str * test_str;
res_unstr = A_unstr * test_unstr;

res_str_interpolated = interp_on_fespace(fespace_str,res_str,fespace_unstr);
err = apply_dirichlet_bc_rhs(abs(res_str_interpolated-res_unstr),fespace_unstr,@(x) [0;0;0;0]);
assert(max(err < 5e-2));


