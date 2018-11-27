function err = compute_H1_error(fespace,values,fexact,gradexact)
% Compute H1 error with respect to exact solution
% input=
%           fespace: finite element space
%           values: vector of degrees of freedom
%           fexact: exact solution
%           gradexact: gradient of the exact solution
%
% output=
%           err = H1 error
%

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;

n_elements = size(connectivity,1);

n_gauss = 6;
[gp,weights,~] = gauss_points2D(n_gauss);

nlocalfunctions = fespace.n_functions_per_element;



err = 0;
if (~strcmp(fespace.mesh.type,'structured'))
    for i = 1:n_elements
        indices = connectivity(i,:);
        x1 = vertices(indices(1),1:2)';
        x2 = vertices(indices(2),1:2)';
        x3 = vertices(indices(3),1:2)';
        
        mattransf = [x2-x1 x3-x1];
        invmat = inv(mattransf);
        
        
        % transformation from parametric to physical
        transf = @(x) mattransf*x + x1;
        dettransf = abs(det(mattransf));
        
        for j = 1:n_gauss
            functions = fespace.functions(gp(:,j));
            grads =  invmat'*fespace.grads(gp(:,j));
            value_in_gp_f = 0;
            value_in_gp_grad = [0;0];
            for k = 1:nlocalfunctions
                value_in_gp_f = value_in_gp_f + values(indices(k))*functions(k);
                value_in_gp_grad = value_in_gp_grad + values(indices(k))*grads(:,k);
            end
            gradex = gradexact(transf(gp(:,j)));
            err = err + dettransf*((value_in_gp_f-fexact(transf(gp(:,j))))^2 + (value_in_gp_grad(1)-gradex(1))^2 + (value_in_gp_grad(2)-gradex(2))^2)*weights(j)/2;
        end
    end
else
    [fespace,gp,weights] = add_members_structured_meshes(fespace,n_gauss);
    for i = 1:n_elements
        indices = connectivity(i,1:end-1);
        x1 = vertices(indices(1),1:2)';
        
        % then the triangle is in this configuration /|
        if (indices(2) == indices(1) + 1)
            for j = 1:n_gauss
                value_in_gp_f = fespace.transffuns1{j}*values(indices);
                value_in_gp_grad = fespace.transfgrads1{j}*values(indices);
                
                gradex = gradexact(fespace.transf1(gp(:,j),x1));
                err = err + fespace.dettransf1*((value_in_gp_f-fexact(fespace.transf1(gp(:,j),x1)))^2 + ...
                    (value_in_gp_grad(1)-gradex(1))^2 + (value_in_gp_grad(2)-gradex(2))^2)*weights(j)/2;
            end
        else
            for j = 1:n_gauss
                value_in_gp_f = fespace.transffuns2{j}*values(indices);
                value_in_gp_grad = fespace.transfgrads2{j}*values(indices);
                
                gradex = gradexact(fespace.transf2(gp(:,j),x1));
                err = err + fespace.dettransf2*((value_in_gp_f-fexact(fespace.transf2(gp(:,j),x1)))^2 + ...
                    (value_in_gp_grad(1)-gradex(1))^2 + (value_in_gp_grad(2)-gradex(2))^2)*weights(j)/2;
            end
        end
    end
end
err = sqrt(err);
