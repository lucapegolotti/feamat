function [T] = assemble_transport(b,fespace)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;
nlocalfunctions = fespace.n_functions_per_element;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);

T = zeros(n_nodes,n_nodes);

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
        transffun = fespace.functions(gp(:,j))';
        transfgrad = invmat' * fespace.grads(gp(:,j));
        for k = 1:nlocalfunctions
            for l = 1:nlocalfunctions
                   transport_element = ((b(transf(gp(:,j)))'*transfgrad(:,k)).*transffun(l))*dettransf*weights(j)/2;
                   T(indices(l),indices(k)) = T(indices(l),indices(k)) + transport_element;
            end
        end
    end
end

T = sparse(T);