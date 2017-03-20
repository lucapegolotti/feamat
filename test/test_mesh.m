clear all
clc

L = 4;
H = 1;

n2 = 5;
n1 = n2*L;

mesh = create_mesh(L,H,n1,n2);
draw_mesh(mesh);
