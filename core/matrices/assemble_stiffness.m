function [A] = assemble_stiffness(mu,fespace)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);

A = sparse(n_nodes,n_nodes);

for i = 1:n_elements
    indices = connectivity(i,1:end-1);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';
    
    mattransf = [x2-x1 x3-x1];
    invmat = inv(mattransf);
    
    % transformation from parametric to physical
    transf = @(x) mattransf*x + x1;       
    dettransf = abs(det(mattransf));
    
    for j = 1:n_gauss
        transfgrad = invmat' * fespace.grads(gp(:,j));
        stiffness_elements = mu(transf(gp(:,j)))*dettransf*(transfgrad'*transfgrad)*weights(j)/2;
        A(indices,indices) = A(indices,indices) + stiffness_elements;
    end
end

A = sparse(A);