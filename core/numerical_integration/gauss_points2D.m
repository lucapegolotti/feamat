function [gp,weights,order] = gauss_points2D(n_gauss)
% Returns gauss quadrature points on a 2D triangle
% input=
%           n_gauss: number of gauss points
%
% output=
%           gp: gauss points
%           weights: weight of each point
%           order: order of the quadrature rule

if (n_gauss == 3)
    gp = [1/6 1/6; 2/3 1/6; 1/6 2/3]';
    weights = [1/3 1/3 1/3];
    order = 2;
elseif (n_gauss == 4)
    gp = [1/3 1/3; 1/5 3/5; 1/5 1/5; 3/5 1/5]';
    weights = [-27/48 25/48 25/48 25/48];
    order = 3;
else
   error(['No 2D Gauss integration rule with ',num2str(n_gauss)',' points!!']); 
end


end

