function [I,code] = evaluate_fe_function(f_dofs,fespace,x_p)
% Evaluate finite element function in point
% ATTENTION! this functions assumes that the underlying mesh is structured
% the function returns code=1 if the point is outside the mesh
% input=
%           f_dofs: values at degrees of freedom of the function          
%           fespace: finite element space
%           x_p: point of interest
% output=
%           I: value of the function
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

mattransf = [x2-x1 x3-x1];

transf = @(xq) mattransf\(xq-x1);    

transfun = fespace.functions(transf([x;y]))';

I = transfun*f_dofs(indv(1:end-1));
code = 0;


