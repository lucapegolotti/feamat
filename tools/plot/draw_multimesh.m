function draw_multimesh(meshes)
% Draw multiple meshes, which must be arranged in a structure
% input=
%           meshes: meshes to draw
%

xp = inf;
yp = inf;
rp = -inf;
up = -inf;

for i = 1:numel(meshes)
    
    mesh = meshes{i};

    n_el = size(mesh.elements,1);

    for j = 1:n_el
        draw_mesh_element(j,mesh.vertices,mesh.elements)
        hold on
    end
    
    xp = min(xp,mesh.xp);
    yp = min(yp,mesh.yp);
    
    rp = max(rp,mesh.xp+mesh.L);
    up = max(up,mesh.yp+mesh.H);
end

axis equal
axis([xp rp yp up]);

end