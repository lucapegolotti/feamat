function err = compute_H1_error(fespace,values,fexact,gradexact)
% Compute L2 error with respect to exact solution
%           fespace: finite element space
%           values: vector of degrees of freedom
%           fexact: exact solution
%

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;

n_elements = size(connectivity,1);

n_gauss = 3;
[gp,weights,~] = gauss_points2D(n_gauss);

nlocalfunctions = fespace.n_functions_per_element;



err = 0;
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
        err = err + dettransf*((value_in_gp_f-fexact(transf(gp(:,j))))^2*0 + (value_in_gp_grad(1)-gradex(1))^2 + (value_in_gp_grad(2)-gradex(2))^2)*weights(j)/2;
    end
end
err = sqrt(err);


