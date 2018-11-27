function [C] = assemble_convective_term(fespace,u_old)
% Assemble connective matrix u_{old} grad u
% input=
%           fespace: finite element space
%           u_old: advecting velocity field
%
% output=
%           C: convective matrix

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);
n_functions = fespace.n_functions_per_element;
n_functionsqr = n_functions^2;

% number of gauss points
n_gauss = 6;

elements_C = zeros(n_functions*n_elements,1);
indices_i = zeros(n_functions*n_elements,1);
indices_j = zeros(n_functions*n_elements,1);

if (~strcmp(fespace.mesh.type,'structured'))
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
            transfun = fespace.functions(gp(:,j));
            
            u1 = transfun'*u_old(indices);
            u2 = transfun'*u_old(indices + n_nodes);
            
            convective_elements = dettransf*transfgrad'*[u1;u2]*transfun'*weights(j)/2;
            new_elements = new_elements + convective_elements(:);
        end
        elements_C(currindices) = new_elements;
    end
else
    [fespace,~,weights] = add_members_structured_meshes(fespace, n_gauss);
    
    for i = 1:n_elements
        indices = connectivity(i,1:end-1);
        
        [I1,I2] = meshgrid(indices,indices);
        currindices = (i-1)*n_functionsqr+1:n_functionsqr*i;
        indices_i(currindices) = I1(:);
        indices_j(currindices) = I2(:);
        
        new_elements = zeros(size(I1,1)^2,1);
        % then the triangle is in this configuration /|
        if (indices(2) == indices(1) + 1)
            for j = 1:n_gauss
                u1 = fespace.transffuns1{j}*u_old(indices);
                u2 = fespace.transffuns1{j}*u_old(indices + n_nodes);
                convective_elements = fespace.dettransf1*fespace.transfgrads1{j}'*[u1;u2]*fespace.transffuns1{j}*weights(j)/2;
                new_elements = new_elements + convective_elements(:);
            end
        else
            for j = 1:n_gauss
                u1 = fespace.transffuns2{j}*u_old(indices);
                u2 = fespace.transffuns2{j}*u_old(indices + n_nodes);
                convective_elements = fespace.dettransf2*fespace.transfgrads2{j}'*[u1;u2]*fespace.transffuns2{j}*weights(j)/2;
                new_elements = new_elements +convective_elements(:);
            end
        end
        elements_C(currindices) = new_elements;
    end
end
C = sparse(indices_i,indices_j,elements_C,n_nodes,n_nodes);
