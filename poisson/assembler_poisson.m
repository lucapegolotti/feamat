function [A,b] = assembler_poisson(fun,mu,fespace,nuemann_functions)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);
bc_flags = fespace.bc;

thereisneumann = 1;

if (length(find(bc_flags)) == 4)
    thereisneumann = 0;
end

A = zeros(n_nodes,n_nodes);
b = zeros(n_nodes,1);

[gp,weights,n_gauss] = gauss_points2D(2);
[gp1d,weights1d,n_gauss1d] = gauss_points1D(3);

nlocalfunctions = fespace.n_functions_per_element;

disp('Building system');

tic 

for i = 1:n_elements
    indices = connectivity(i,:);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';
    
    % if boundary element, find vertices that are on edges
    if (thereisneumann)
        if (connectivity(i,end) > 0)
            if (vertices(indices(1),3) > 0 && vertices(indices(2),3) > 0)
                x1_edge = x1;
                x2_edge = x2;
                transfgp_1d_to_2d = @(x) [1;0]/2*x + [1;0]/2;
            elseif (vertices(indices(2),3) > 0 && vertices(indices(3),3) > 0)
                x1_edge = x2;
                x2_edge = x3;
                transfgp_1d_to_2d = @(x) [-1;1]/2*x + [1;1]/2;
            elseif (vertices(indices(3),3) > 0 && vertices(indices(1),3) > 0)
                x1_edge = x3;
                x2_edge = x1;
                transfgp_1d_to_2d = @(x) [0;-1]/2*x + [0;1]/2;
            else
                error('An error occured in boundary element detection!');
            end
        end
    end

    mattransf = [x2-x1 x3-x1];
    invmat = inv(mattransf);
    
    % transformation from parametric to physical
    transf = @(x) mattransf*x + x1;       
    dettransf = abs(det(mattransf));
    
    for j = 1:n_gauss
        transfgrad = invmat' * fespace.grads(gp(:,j));
        functions = fespace.functions(gp(:,j));
        for k = 1:nlocalfunctions
            b(indices(k)) = b(indices(k)) + dettransf*fun(transf(gp(:,j)))*functions(k)*weights(j)/2;
            for l = 1:nlocalfunctions
                stiffness_element = mu(transf(gp(:,j)))*dettransf*transfgrad(:,k)'*transfgrad(:,l)*weights(j)/2;
                A(indices(k),indices(l)) = A(indices(k),indices(l)) + stiffness_element;
            end
        end
    end
    
    if (thereisneumann)
        transf1d = @(x) (x2_edge - x1_edge)/2*x + (x1_edge + x2_edge)/2;
        dtransf = norm((x2_edge - x1_edge)/2);
        for j = 1:n_gauss1d
            functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
            for k = 1:nlocalfunctions
                if (nodes(indices(k),3)>0 && bc_flags(nodes(indices(k),3))==0)
                    nf = nuemann_functions(transf1d(gp1d(j)));
                    b(indices(k)) = b(indices(k)) - dtransf*nf(nodes(indices(k),3))*functions(k)*weights1d(j);
                end
            end
        end
    end
end


elapsed = toc;
disp(['Elapsed time = ', num2str(elapsed),' s']);
disp('------------------------------');


