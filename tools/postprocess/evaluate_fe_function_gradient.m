function [I,code] = evaluate_fe_function_gradient(f_dofs,fespace,x_p)
% Evaluate gradient of finite element function in point
% input=
%           f_dofs: values at degrees of freedom of the function          
%           fespace: finite element space
%           x_p: point of interest
% output=
%           I: value of the gradient of the function
%           code: 0 if ok, 1 if the point is outside the domain
mesh = fespace.mesh;

x = x_p(1);
y = x_p(2);
if (x < mesh.xp || x > mesh.xp+mesh.L)
    I = 0;
    code = 1;
    return;
end
if (y < mesh.yp || y > mesh.yp+mesh.H)
    I = 0;
    code = 1;
    return;
end

if (strcmp(mesh.type,'structured'))

    % find element that contains (x,y)
    xx = mesh.X(1:end-1,1);
    yy = mesh.Y(1,1:end-1);

    L = length(xx);

    indx = max(find(xx <= x));
    indy = max(find(yy <= y));

    elindex = 2*L*(indy-1) + indx;

    % check if (x,y) is in this element or the one above
    indv = fespace.connectivity(elindex,:);

    x1 = mesh.vertices(indv(1),1:2)';
    x2 = mesh.vertices(indv(2),1:2)';
    x3 = mesh.vertices(indv(3),1:2)';

    m1 = (x3(2)-x1(2))/(x3(1)-x1(1));
    m2 = (y-x1(2))/(x-x1(1));

    % then we chose the wrong element: choosing the one above
    if (m1 < m2)
        indv = fespace.connectivity(elindex+L,:);
        x1 = mesh.vertices(indv(1),1:2)';
        x2 = mesh.vertices(indv(2),1:2)';
        x3 = mesh.vertices(indv(3),1:2)';
    end
elseif (strcmp(mesh.type,'unstructured'))
    aux = mesh.vertices(:,1:2) - x_p(:)';
    distances = sqrt(aux(:,1).^2 + aux(:,2).^2);
    
    % search closest vertex to input point
    [~,index] = min(distances);

    elements_with_vertex = (mesh.elements(:,1) == index) + ...
                           (mesh.elements(:,2) == index) + ...
                           (mesh.elements(:,3) == index);
    indices_elements = find(elements_with_vertex);
    
    found = 0;
    for i = 1:size(indices_elements,1)
        idx = indices_elements(i);
        x1 = mesh.vertices(mesh.elements(idx,1),1:2)';
        x2 = mesh.vertices(mesh.elements(idx,2),1:2)';
        x3 = mesh.vertices(mesh.elements(idx,3),1:2)';

        P12 = (x1-x2)'; P23 = (x2-x3)'; P31 = (x3-x1)';
        
        s1 = sign(det([P31;P23]))*sign(det([x3'-x_p(:)';P23]));
        s2 = sign(det([P12;P31]))*sign(det([x1'-x_p(:)';P31]));
        s3 = sign(det([P23;P12]))*sign(det([x2'-x_p(:)';P12]));
        is_inside =  s1 >= -1e-16 * s2 >= -1e-16 * s3 > -1e-16;
        if (is_inside)
            found = 1;
            break;
        end
    end
    if (found == 0)
        I = 0;
        code = 1;
        return;
    end
    indv = fespace.connectivity(idx,:);
end

mattransf = [x2-x1 x3-x1];
invmat = inv(mattransf);

transf = @(xq) mattransf\(xq-x1);    

transfgrad = fespace.grads(transf([x;y]));

I = invmat'*transfgrad*f_dofs(indv(1:end-1));
code = 0;


