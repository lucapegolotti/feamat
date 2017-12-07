% common variables

tol = 1e-12;

%% Test 1: all fields assigned correctly
xp = 0;
yp = 0;

L = 1;
H = 1;

mesh = create_mesh(xp,yp,L,H,3,3);

assert(mesh.xp == 0);
assert(mesh.yp == 0); 
assert(mesh.L == L);
assert(mesh.H == 1);
assert(strcmp(mesh.type,'structured'));
assert(mesh.h == 1/3)

%% Test 2: verify vertex positions

xp = 0;
yp = 0;
L = 2;
H = 1;

mesh = create_mesh(xp,yp,L,H,2,2);

assert(norm(mesh.vertices(1,:)-[0 0 1 4]) < tol);
assert(norm(mesh.vertices(2,:)-[1 0 1 0]) < tol);
assert(norm(mesh.vertices(5,:)-[1 0.5 0 0]) < tol);
assert(norm(mesh.vertices(end,:)-[2 1 2 3]) < tol);

%% Test 3: verify elements number and ordering
xp = 0;
yp = 0;
L = 1;
H = 1;

mesh = create_mesh(xp,yp,L,H,4,4);

assert(size(mesh.elements,1) == 4*4*2);
assert(size(mesh.elements,2) == 4);
assert(norm(mesh.elements(1,:)-[1 2 7 1]) < tol);

%% Test 4: test mesh with one element
xp = 0.5;
yp = 0.5;
L = 0.5;
H = 1;

mesh = create_mesh(xp,yp,L,H,1,1);

assert(size(mesh.elements,1) == 1*1*2)

