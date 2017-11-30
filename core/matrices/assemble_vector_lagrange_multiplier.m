function [v,vt] = assemble_vector_lagrange_multiplier(fespace)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);

v = zeros(n_nodes,1);

for i = 1:n_elements
    indices = connectivity(i,1:end-1);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';
    
    mattransf = [x2-x1 x3-x1];
    
    % transformation from parametric to physical
    dettransf = abs(det(mattransf));
    
    for j = 1:n_gauss
        transffun = fespace.functions(gp(:,j))';
        element = dettransf*(transffun)*weights(j)/2;
        v(indices) = v(indices) + element';
    end
end

vt = v';
v = apply_dirichlet_bc_rhs(v,fespace,@(x) [0;0;0;0]);
