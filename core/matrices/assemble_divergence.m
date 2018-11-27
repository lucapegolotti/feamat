function [B] = assemble_divergence(fespace_u,fespace_p,derivative)

if (derivative == 'dx')
    index_der = 1;
elseif (derivative == 'dy')
    index_der = 2;
else
    error([derivative, ' is not a valid derivative!']);
end

connectivity_u = fespace_u.connectivity;
connectivity_p = fespace_p.connectivity;

vertices = fespace_u.mesh.vertices;
nodes_u = fespace_u.nodes;
nlocalfunctions_u = fespace_u.n_functions_per_element;

nodes_p = fespace_p.nodes;
nlocalfunctions_p = fespace_p.n_functions_per_element;

n_functionsqr = nlocalfunctions_u*nlocalfunctions_p;

n_elements = size(connectivity_u,1);

n_nodes_u = size(nodes_u,1);
n_nodes_p = size(nodes_p,1);

n_gauss = 6;
[gp,weights,~] = gauss_points2D(n_gauss);

elements_B = zeros(nlocalfunctions_u*nlocalfunctions_p*n_elements,1);
indices_i = zeros(nlocalfunctions_u*nlocalfunctions_p*n_elements,1);
indices_j = zeros(nlocalfunctions_u*nlocalfunctions_p*n_elements,1);


if (~strcmp(fespace_u.mesh.type,'structured'))
    for i = 1:n_elements
        indices_u = connectivity_u(i,1:end-1);
        indices_p = connectivity_p(i,1:end-1);        
        x1 = vertices(indices_u(1),1:2)';
        x2 = vertices(indices_u(2),1:2)';
        x3 = vertices(indices_u(3),1:2)';

        [I1,I2] = meshgrid(indices_u,indices_p);
        currindices = (i-1)*n_functionsqr+1:n_functionsqr*i;
        indices_i(currindices) = I2(:);
        indices_j(currindices) = I1(:);
        
        mattransf = [x2-x1 x3-x1];
        invmat = inv(mattransf);
        
        % transformation from parametric to physical
        dettransf = abs(det(mattransf));

        new_elements = zeros(n_functionsqr,1);
        
        for j = 1:n_gauss
            transfgrad = invmat'*fespace_u.grads(gp(:,j));
            transfun = fespace_p.functions(gp(:,j));

            divergence_element = dettransf*transfun(:)* ...
                                 (transfgrad(index_der,:))*weights(j)/2;
            new_elements = new_elements + divergence_element(:);
        end
        elements_B(currindices) = new_elements;
    end
else
    [fespace_u,~] = add_members_structured_meshes(fespace_u, n_gauss);
    [fespace_p,~] = add_members_structured_meshes(fespace_p, n_gauss);
    
    divergence_elements1 = {};
    divergence_elements2 = {};
    for j = 1:n_gauss
        divergence_elements1{end+1} = fespace_u.dettransf1* ...
            fespace_p.transffuns1{j}(:)* ...
            (fespace_u.transfgrads1{j}(index_der,:))*weights(j)/2;

        divergence_elements2{end+1} = fespace_u.dettransf2* ...
            fespace_p.transffuns2{j}(:)* ...
            (fespace_u.transfgrads2{j}(index_der,:))*weights(j)/2;
    end
    sum_divergence_elements1 = divergence_elements1{1};
    sum_divergence_elements2 = divergence_elements2{1};
    for i = 2:n_gauss
        sum_divergence_elements1 = sum_divergence_elements1 + divergence_elements1{i};
        sum_divergence_elements2 = sum_divergence_elements2 + divergence_elements2{i};
    end
    
    for i = 1:n_elements
        indices_u = connectivity_u(i,1:end-1);
        indices_p = connectivity_p(i,1:end-1);
        
        [I1,I2] = meshgrid(indices_u,indices_p);
        currindices = (i-1)*n_functionsqr+1:n_functionsqr*i;
        indices_i(currindices) = I2(:);
        indices_j(currindices) = I1(:);
        
        % then the triangle is in this configuration /|
        if (indices_u(2) == indices_u(1) + 1)
            new_elements = sum_divergence_elements1(:);
        else
            new_elements = sum_divergence_elements2(:);
        end
        elements_B(currindices) = new_elements;
    end
end
B = sparse(indices_i,indices_j,elements_B,n_nodes_p,n_nodes_u);
