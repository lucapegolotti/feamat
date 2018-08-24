function [mesh,fespace] = move_mesh(fespace,vec)
% Applies the prescribed displacement field to all vertices and nodes of
% the mesh
% input=
%           fespace: finite elemnet space of displacement
%           vec: dofs of the displacement
% output=
%           mesh: modified mesh
%           fespace: fespace with modified mesh

% retrieve mesh from fespace
mesh = fespace.mesh;

n_vertices = size(mesh.vertices,1);

u1 = vec.u1;
u2 = vec.u2;

mesh.vertices(:,1) = mesh.vertices(:,1) + u1(1:n_vertices);
mesh.vertices(:,2) = mesh.vertices(:,2) + u2(1:n_vertices);

fespace.nodes(:,1) = fespace.nodes(:,1) + u1;
fespace.nodes(:,2) = fespace.nodes(:,2) + u2;

% update fields in mesh
mesh.xp = min(mesh.vertices(:,1));
mesh.yp = min(mesh.vertices(:,2));
mesh.L = max(mesh.vertices(:,1)) - mesh.xp;
mesh.H = max(mesh.vertices(:,2)) - mesh.yp;

n_elements = size(mesh.elements,1);

h = 0;
% update max h
for i = 1:n_elements
    i1 = mesh.elements(i,1);
    i2 = mesh.elements(i,2);
    i3 = mesh.elements(i,3);
    
    v1 = mesh.vertices(i1,1:2)';
    v2 = mesh.vertices(i2,1:2)';
    v3 = mesh.vertices(i3,1:2)';
    
    h1 = norm(v1 - v2);
    h2 = norm(v2 - v3);
    h3 = norm(v3 - v1);
    
    h = max(h,h1);
    h = max(h,h2);
    h = max(h,h3);
end


fespace.mesh = mesh;
