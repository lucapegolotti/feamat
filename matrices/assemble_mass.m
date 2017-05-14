function [M] = assemble_mass(fespace)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;
nlocalfunctions = fespace.n_functions_per_element;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);

M = zeros(n_nodes,n_nodes);

for i = 1:n_elements
    indices = connectivity(i,:);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';
    
    mattransf = [x2-x1 x3-x1];
    
    % transformation from parametric to physical
    dettransf = abs(det(mattransf));
    
    for j = 1:n_gauss
        transffun = fespace.functions(gp(:,j))';
        for k = 1:nlocalfunctions
            for l = 1:nlocalfunctions
                mass_element = dettransf*transffun(k)'*transffun(l)*weights(j)/2;
                M(indices(k),indices(l)) = M(indices(k),indices(l)) + mass_element;
            end
        end
    end
end

M = sparse(M);