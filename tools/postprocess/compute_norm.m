function norm = compute_norm(fespace,values,type)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;

n_elements = size(connectivity,1);

n_gauss = 3;
[gp,weights,order] = gauss_points2D(n_gauss);

nlocalfunctions = fespace.n_functions_per_element;


if (type == 'L2')
    norm = 0;
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
                norm = norm + dettransf*values(indices(k))^2*functions(k)^2*weights(j)/2;
            end
        end
    end  
    norm = sqrt(norm);

elseif (type == 'H1')
    
    norm = 0;
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
            grads = fespace.grads(gp(:,j));
            for k = 1:nlocalfunctions
                norm = norm + dettransf*values(indices(k))^2*(functions(k)^2 + grads(1,k)^2 + grads(2,k)^2)*weights(j)/2;
            end
        end
    end  
    norm = sqrt(norm);
else
    error([type,' norm is not implemented!']);
end


end

