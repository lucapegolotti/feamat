function [A,b] = assembler_poisson(fun,mu,fespace,dirichlet_functions)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

A = zeros(n_nodes,n_nodes);

gp = [1/6 1/6; 2/3 1/6; 1/6 2/3]';
weights = [1/3 1/3 1/3];
n_gauss = length(weights);

nlocalfunctions = size(fespace.grads,2);

mu_vec = @(x) mu(x(1),x(2));

disp('Building system');

for i = 1:n_elements
    indices = connectivity(i,1:3);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';

    transf = [x2-x1 x3-x1];
    invtransf = inv(transf);
    dettransf = abs(det(transf));
    
    transfgrad = invtransf' * fespace.grads;

    for j = 1:n_gauss
        for k = 1:nlocalfunctions
            for l = 1:nlocalfunctions
                A(indices(k),indices(l)) = A(indices(k),indices(l)) + ...
                mu_vec(invtransf*gp(:,j))*dettransf*transfgrad(:,k)'*transfgrad(:,l)*weights(j);
            end
        end
    end
end

b = zeros(n_nodes,1);

for i = 1:n_nodes
    b(i) = fun(nodes(i,1),nodes(i,2));
end

[A,b] = apply_bc(A,b,nodes,fespace.bc,dirichlet_functions);


