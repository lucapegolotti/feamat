function [A] = assemble_stiffness(mu,fespace)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;
nlocalfunctions = fespace.n_functions_per_element;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);

A = zeros(n_nodes,n_nodes);

for i = 1:n_elements
    indices = connectivity(i,:);
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
        for k = 1:nlocalfunctions
            for l = 1:nlocalfunctions
                stiffness_element = mu(transf(gp(:,j)))*dettransf*transfgrad(:,k)'*transfgrad(:,l)*weights(j)/2;
                A(indices(k),indices(l)) = A(indices(k),indices(l)) + stiffness_element;
            end
        end
    end
end

A = sparse(A);