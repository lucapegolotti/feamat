function [A,b] = assembler_poisson(vertices, connectivity,flags_bc)

n = size(vertices,1);

A = zeros(n,n);

gp = [1/6 1/6; 2/3 1/6; 1/6 2/3];
weights = [1/3 1/3 1/3];
n_gauss = length(weights);

grad1_ref = [-1;-1];
grad2_ref = [1;0];
grad3_ref = [0;1];

for i = 1:n
    indices = connectivity(i,1:3);
    x1 = vertices(indices(1),:)';
    x2 = vertices(indices(2),:)';
    x3 = vertices(indices(3),:)';

    transf = [x2-x1 x3-x1];
    invtransf = inv(transf);
    dettransf = abs(det(transf));
    
    grad1 = invtransf' * grad1_ref
    grad2 = invtransf' * grad2_ref
    grad3 = invtransf' * grad3_ref

    for j = 1:n_gauss
        A(indices(1),indices(2)) = A(indices(1),indices(2)) + dettransf * grad1'*grad2*weights(j) / 2;
        A(indices(2),indices(3)) = A(indices(2),indices(3)) + dettransf * grad2'*grad3*weights(j) / 2;
        A(indices(3),indices(1)) = A(indices(3),indices(1)) + dettransf * grad3'*grad1*weights(j) / 2;
        A(indices(1),indices(1)) = A(indices(1),indices(1)) + dettransf * grad1'*grad1*weights(j) / 2;
        A(indices(2),indices(2)) = A(indices(2),indices(2)) + dettransf * grad2'*grad2*weights(j) / 2;
        A(indices(3),indices(3)) = A(indices(3),indices(3)) + dettransf * grad3'*grad3*weights(j) / 2;
    end
    A(indices(2),indices(1)) = A(indices(1),indices(2));
    A(indices(3),indices(2)) = A(indices(2),indices(3));
    A(indices(1),indices(3)) = A(indices(3),indices(1));
end


