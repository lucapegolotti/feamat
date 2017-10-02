clear all
clc

% create mesh on left subdomain
xp1 = 0;
yp1 = 0;
L1 = 0.5;
H1 = 1;

mesh1 = create_mesh(xp1,yp1,L1,H1,10,20);

% create mesh on right subdomain
xp2 = 0.5;
yp2 = 0;
L2 = 0.5;
H2 = 1;

mesh2 = create_mesh(xp2,yp2,L2,H2,20,20);

% draw global mesh
meshes = {};
meshes{end+1} = mesh1;
meshes{end+1} = mesh2;
draw_multimesh(meshes);
