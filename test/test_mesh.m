%% test creation of single mesh
clear all
clc

L = 4;
H = 1;

n2 = 5;
n1 = n2*L;

mesh = create_mesh(0,0,L,H,n1,n2);
draw_mesh(mesh);

%% test creation of multimesh

clear all
clc

L = 4;
H = 1;

n2 = 5;
n1 = n2*L;

mesh1 = create_mesh(0,0,L,H,n1,n2);


L = 2;
H = 1;

n2 = 5;
n1 = n2*L*2;

mesh2 = create_mesh(4,0,L,H,n1,n2);

meshes = {};
meshes{end+1} = mesh1;
meshes{end+1} = mesh2;

draw_multimesh(meshes);




