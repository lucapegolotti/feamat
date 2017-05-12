function [b] = assemble_rhs(fespace,fun)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;
nlocalfunctions = fespace.n_functions_per_element;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);

b = zeros(n_nodes,1);

for i = 1:n_elements
    indices = connectivity(i,:);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';
    
    mattransf = [x2-x1 x3-x1];
    
    % transformation from parametric to physical
    transf = @(x) mattransf*x + x1;       
    dettransf = abs(det(mattransf));
    
    for j = 1:n_gauss
        functions = fespace.functions(gp(:,j));
        for k = 1:nlocalfunctions
            b(indices(k)) = b(indices(k)) + dettransf*fun(transf(gp(:,j)))*functions(k)*weights(j)/2;
        end
    end
end