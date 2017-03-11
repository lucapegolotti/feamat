clear all
clc

L = 4;
H = 1;

n2 = 20;
n1 = n2*L;

[X,Y,vert,conn] = create_mesh(L,H,n1,n2);

n_el = size(conn,1);

for i = 1:n_el
    draw_mesh_element(i,vert,conn)
    hold on
end
axis equal
axis([0 L 0 H]);
