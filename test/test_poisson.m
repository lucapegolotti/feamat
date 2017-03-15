clear all
close all
clc

L = 2;
H = 1;

n2 = 30;
n1 = n2*L;

[X,Y,vert,conn] = create_mesh(L,H,n1,n2);

n_el = size(conn,1);

for i = 1:n_el
    draw_mesh_element(i,vert,conn)
    hold on
end
axis equal
axis([0 L 0 H]);

f = @(x,y) 1;

[A,b] = assembler_poisson(f,vert,conn,[1 1 1 1]);

figure()
spy(A)

sol = A\b;

n1 = size(X,1);
n2 = size(X,2);


surf(X,Y,reshape(sol,n1,n2));

% mat = [-3 3; 3 0];
% grad = mat*[0;1];
% hold on
% plot([0 grad(2)]+1/3,[0 grad(1)]+1/3)


