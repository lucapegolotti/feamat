function [B] = assemble_divergence(fespace_u,fespace_p,derivative)

if (derivative == 'dx')
    index_der = 1;
elseif (derivative == 'dy')
    index_der = 2;
else
   error([derivative, ' is not a valid derivative!']); 
end

connectivity_u = fespace_u.connectivity;
vertices = fespace_u.mesh.vertices;
nodes_u = fespace_u.nodes;
nlocalfunctions_u = fespace_u.n_functions_per_element;

nodes_p = fespace_p.nodes;
nlocalfunctions_p = fespace_p.n_functions_per_element;

n_elements = size(connectivity_u,1);

n_nodes_u = size(nodes_u,1);
n_nodes_p = size(nodes_p,1);

[gp,weights,n_gauss] = gauss_points2D(2);

B = sparse(n_nodes_p,n_nodes_u);

for i = 1:n_elements
    indices = connectivity_u(i,1:end-1);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';
    
    mattransf = [x2-x1 x3-x1];
    invmat = inv(mattransf);
    
    % transformation from parametric to physical
    dettransf = abs(det(mattransf));
    
    for j = 1:n_gauss
        transfgrad = invmat'*fespace_u.grads(gp(:,j));
        transfun = fespace_p.functions(gp(:,j));
        for k = 1:nlocalfunctions_u
            for l = 1:nlocalfunctions_p
                divergence_element = dettransf*(transfgrad(index_der,k))*transfun(l)*weights(j)/2;
                B(indices(l),indices(k)) = B(indices(l),indices(k)) + divergence_element;
            end
        end
    end
end
