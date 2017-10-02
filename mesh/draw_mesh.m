function draw_mesh(mesh)

n_el = size(mesh.elements,1);

for i = 1:n_el
    draw_mesh_element(i,mesh.vertices,mesh.elements)
    hold on
end

axis equal
axis([mesh.xp mesh.xp+mesh.L mesh.yp mesh.yp+mesh.H]);

end