function [I,code] = evaluate_fe_function_gradient(f_dofs,fespace,x_p)
% Evaluate gradient of finite element function in point or in points
% input=
%           f_dofs: values at degrees of freedom of the function
%           fespace: finite element space
%           x_p: point(s) of interest. If x_p is a matrix of points,
%                the size must be 2 x (nb points).
% output=
%           I: value(s) of the gradient of the function.
%              Matrix of size 2 x nb points.
%           code: 0 if ok, 1 if the point is outside the domain.
%                 If multiple points, vector of size nb points.

nb_points = size(x_p,2);
I = zeros(2,nb_points);
code = zeros(nb_points,1);
for i = 1:size(x_p,2)
    
    [index,x1,x2,x3] = find_element_containing_point(fespace.mesh,x_p(:,i));
    
    % perturb the point as it may lay on a line
    if (index == -1)
        count = 0;
        rot = [0 -1; 1 0];
        pert = [1e-8;0];
        while (index == -1 && count < 4)
            count = count + 1;
            x_p_new = x_p(:,i) + pert;
            [index,x1,x2,x3] = find_element_containing_point(fespace.mesh,x_p_new);
            pert = rot*pert;
        end
        
        pert = [1e-8;1e-8];
        count = 0;
        while (index == -1 && count < 4)
            count = count + 1;
            x_p_new = x_p(:,i) + pert;
            [index,x1,x2,x3] = find_element_containing_point(fespace.mesh,x_p_new);
            pert = rot*pert;
        end
    end
    
    % perturb the point as it may lay on a line
    count = 0;
    while (index == -1 && count < 4)
        count = count + 1;
        if (count < 3)
            x_p_new = x_p(:,i) + [(-1)^count * 1e-5;0];
        else
            x_p_new = x_p(:,i) + [0;(-1)^count * 1e-5];
        end
        [index,x1,x2,x3] = find_element_containing_point(fespace.mesh,x_p_new);
    end
    
    indv = fespace.connectivity(index,:);
    
    mattransf = [x2-x1 x3-x1];
    invmat = inv(mattransf);
    
    transf = @(xq) mattransf\(xq-x1);
    
    transfgrad = fespace.grads(transf(x_p(:,i)));
    
    I(:,i) = (invmat'*transfgrad*f_dofs(indv(1:end-1)))';
    code(i) = 0;
end


