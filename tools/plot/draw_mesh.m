function draw_mesh(mesh,varargin)
% Draw the mesh (warning: this function is yet to be optimized)
% input=
%           mesh: mesh to draw
%           (optional)
%           'color' to color the boundary elements with different colors
%

color = 0;
if (nargin == 2)
    if (strcmp(varargin{1},'color'))
        color = 1;
    else
        error('Unknown parameter to draw_mesh function!');
    end
end

n_el = size(mesh.elements,1);

if (~color)
    for i = 1:n_el
        draw_mesh_element(i,mesh.vertices,mesh.elements)
        hold on
    end
elseif (color)
    for i = 1:n_el
        draw_mesh_element(i,mesh.vertices,mesh.elements,'color')
        hold on
    end
end

axis equal
axis([mesh.xp mesh.xp+mesh.L mesh.yp mesh.yp+mesh.H]);

end