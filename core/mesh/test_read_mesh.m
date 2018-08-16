% common variables

tol = 1e-2;

%% Test 1: all fields assigned correctly
xp = 0;
yp = 0;

boundary_indicators = @(x) [(norm(x) < 0.3 + tol);
                            abs(x(1)-0.5) < tol;
                            abs(x(1)+0.5) < tol;
                            abs(x(2)-0.5) < tol;
                            abs(x(2)+0.5) < tol];

mesh = read_mesh('../../data/square_hole.msh',boundary_indicators);

assert(mesh.xp == -0.5);
assert(mesh.yp == -0.5); 
assert(mesh.L == 1);
assert(mesh.H == 1);
assert(abs(mesh.h-0.1086) < tol)