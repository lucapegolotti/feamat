function [C] = assemble_convective_term(fespace,u_old)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;
nlocalfunctions = fespace.n_functions_per_element;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);

C = zeros(n_nodes,n_nodes);

for i = 1:n_elements
    indices = connectivity(i,:);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';
    
    mattransf = [x2-x1 x3-x1];
    invmat = inv(mattransf);
    
    % transformation from parametric to physical
    dettransf = abs(det(mattransf));
    
    for j = 1:n_gauss
        transfgrad = invmat' * fespace.grads(gp(:,j));
        transfun = fespace.functions(gp(:,j));
        u1 = 0;
        u2 = 0;
        for k = 1:nlocalfunctions
            u1 = u1 + transfun(k)*u_old(indices(k));
            u2 = u2 + transfun(k)*u_old(indices(k)+n_nodes);
        end
        for k = 1:nlocalfunctions
            for l = 1:nlocalfunctions
                convective_element = dettransf*transfgrad(:,k)'*[u1;u2]*weights(j)/2;
                C(indices(k),indices(l)) = C(indices(k),indices(l)) + convective_element;
            end
        end
    end
end

C = sparse(C);