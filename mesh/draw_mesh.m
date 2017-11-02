function draw_mesh(mesh,varargin)

n_el = size(mesh.elements,1);

if (nargin == 1)
    for i = 1:n_el
        draw_mesh_element(i,mesh.vertices,mesh.elements)
        hold on
    end
elseif (strcmp(varargin{1},'nocolor'))
    for i = 1:n_el
        draw_mesh_element(i,mesh.vertices,mesh.elements,varargin{1})
        hold on
    end
else
    error('Optional parameter is not supported!');
end

axis equal
axis([mesh.xp mesh.xp+mesh.L mesh.yp mesh.yp+mesh.H]);

end