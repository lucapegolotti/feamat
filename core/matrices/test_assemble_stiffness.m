% common variables

tol = 1e-12;

%% Test 1: stiffness on few elements

mesh = create_mesh(0,0,1,1,50,50);
fespace = create_fespace(mesh,'P2',[0 0 0 0]);

tic
A1 = assemble_stiffness(@(x) 1, fespace,1);
toc

fespace.mesh.type = '';
tic
A2 = assemble_stiffness(@(x) 1, fespace);
toc
