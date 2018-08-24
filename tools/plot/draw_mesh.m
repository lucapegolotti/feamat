function draw_mesh(mesh,varargin)
% Draw the mesh (warning: this function is yet to be optimized)
% input=
%           mesh: mesh to draw
%

if (nargin == 1)
    color = [0.3 0.3 0.3];
else
    color = varargin{1};
end

trimesh(mesh.elements,mesh.vertices(:,1),mesh.vertices(:,2),'Color',color);

axis equal
axis([mesh.xp mesh.xp+mesh.L mesh.yp mesh.yp+mesh.H]);