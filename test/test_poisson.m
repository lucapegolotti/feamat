clear all
close all
clc

L = 1;
H = 1;

n2 = 8;
n1 = n2*L;

[X,Y,vert,conn] = create_mesh(L,H,n1,n2);

n_el = size(conn,1);

for i = 1:n_el
    draw_mesh_element(i,vert,conn)
    hold on
end
axis equal
axis([0 L 0 H]);

A = assembler_poisson(vert,conn,1);

figure()
spy(A)

% mat = [-3 3; 3 0];
% grad = mat*[0;1];
% hold on
% plot([0 grad(2)]+1/3,[0 grad(1)]+1/3)


