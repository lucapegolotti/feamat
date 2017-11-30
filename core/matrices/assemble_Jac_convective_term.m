function [C] = assemble_Jac_convective_term(fespace,u_old,blocki,blockj)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;
nlocalfunctions = fespace.n_functions_per_element;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);

C = zeros(n_nodes,n_nodes);

for i = 1:n_elements
    indices = connectivity(i,1:end-1);
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
        du1 = 0;
        du2 = 0;
        for k = 1:nlocalfunctions
            du1 = du1 + transfgrad(:,k)*u_old(indices(k));
            du2 = du2 + transfgrad(:,k)*u_old(indices(k)+n_nodes);
        end
        du = [du1 du2]';
        
        for k = 1:nlocalfunctions
            for l = 1:nlocalfunctions
                mult = dettransf*weights(j)/2;
                convective_element_jac = mult*(du(blocki,blockj))*transfun(l)*transfun(k);
                C(indices(l),indices(k)) = C(indices(l),indices(k)) + convective_element_jac;
            end
        end
    end
end

C = sparse(C);