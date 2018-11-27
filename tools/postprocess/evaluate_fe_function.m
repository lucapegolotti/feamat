function [I,code] = evaluate_fe_function(f_dofs,fespace,x_p)
% Evaluate finite element function in point or multiple points
% input=
%           f_dofs: values at degrees of freedom of the function
%           fespace: finite element space
%           x_p: point(s) of interest. If x_p is a matrix of points,
%                the size must be 2 x (nb points).
% output=
%           I: value(s) of the function. Vector of size nb points.
%           code: 0 if ok, 1 if the point is outside the domain.
%                 If multiple points, vector of size nb points.
if (size(x_p,1) ~= 2)
    x_p = x_p';
end
nb_points = size(x_p,2);
I = zeros(nb_points,1);
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
    
    if (index == -1)
        I = 0;
        code = 1;
        return
    end
    
    indv = fespace.connectivity(index,:);
    mattransf = [x2-x1 x3-x1];
    
    transf = @(xq) mattransf\(xq-x1);
    transfun = fespace.functions(transf(x_p(:,i)))';
    
    I(i) = transfun*f_dofs(indv(1:end-1));
    code(i) = 0;
end


