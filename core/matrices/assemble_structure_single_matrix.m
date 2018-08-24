function [A] = assemble_structure_single_matrix(fespace,block_i,block_j)
% Assemble a single block in a linear elasticity problem
% input=
%           fespace: finite element space
%           block_i: index of the first component
%           block_j: index of the second component
% output=
%           A: output matrix (sparse)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);
n_functions = fespace.n_functions_per_element;
n_functionsqr = n_functions^2;

% number of gauss points
n_gauss = 3;

elements_A = zeros(n_functions*n_elements,1);
indices_i = zeros(n_functions*n_elements,1);
indices_j = zeros(n_functions*n_elements,1);
    
[gp,weights,~] = gauss_points2D(n_gauss);

for i = 1:n_elements
    indices = connectivity(i,1:end-1);
    x1 = vertices(indices(1),1:2)';
    x2 = vertices(indices(2),1:2)';
    x3 = vertices(indices(3),1:2)';

    [I1,I2] = meshgrid(indices,indices);

    currindices = (i-1)*n_functionsqr+1:n_functionsqr*i;
    indices_i(currindices) = I1(:);
    indices_j(currindices) = I2(:);

    new_elements = zeros(size(I1,1)^2,1);

    mattransf = [x2-x1 x3-x1];
    invmat = inv(mattransf);

    % transformation from parametric to physical
    dettransf = abs(det(mattransf));

    for j = 1:n_gauss
        transfgrad = invmat' * fespace.grads(gp(:,j));
        transfgrad_i = transfgrad(block_i,:);
        transfgrad_j = transfgrad(block_j,:);
        stiffness_elements = dettransf*(transfgrad_j'*transfgrad_i)*weights(j)/2;
        new_elements = new_elements + stiffness_elements(:);
    end
    elements_A(currindices) = new_elements;
end

A = sparse(indices_i,indices_j,elements_A,n_nodes,n_nodes);


