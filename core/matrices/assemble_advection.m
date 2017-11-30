function [M] = assemble_advection(c,fespace)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;
nlocalfunctions = fespace.n_functions_per_element;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);

M = zeros(n_nodes,n_nodes);

for i = 1:n_elements
    indices = connectivity(i,1:end-1);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';
    
    mattransf = [x2-x1 x3-x1];
    
    % transformation from parametric to physical
    transf = @(x) mattransf*x + x1;       
    dettransf = abs(det(mattransf));
    
    for j = 1:n_gauss
        transffun = fespace.functions(gp(:,j))';
        advection_element = c(transf(gp(:,j)))*dettransf*(transffun'*transffun)*weights(j)/2;
        M(indices,indices) = M(indices,indices) + advection_element;
    end
end

M = sparse(M);