function [gp,weights,order] = gauss_points1D(n_gauss)
% Returs gauss quadrature points in 1D on the segmente (-1,1)
% input=
%           n_gauss: number of gauss points
%
% output=
%           gp: gauss points
%           weights: weight of each point
%           order: order of the quadrature rule

if (n_gauss == 2)
    gp = [-sqrt(1/3); sqrt(1/3)]';
    weights = [1 1];
    order = 3;
elseif (n_gauss == 3)
    gp = [-sqrt(3/5);0;sqrt(3/5)]';
    weights = [5/9 8/9 5/9];
    order = 5;
elseif (n_gauss == 4)
    gp = [-sqrt(3/7-2/7*sqrt(6/5));+sqrt(3/7-2/7*sqrt(6/5));-sqrt(3/7+2/7*sqrt(6/5));+sqrt(3/7+2/7*sqrt(6/5))]';
    weights = [18+sqrt(30) 18+sqrt(30) 18-sqrt(30) 18-sqrt(30)]/36;
    order = 7;
else
   error(['No 1D Gauss integration rule with',num2str(n_gauss)',' points!!']); 
end


end



