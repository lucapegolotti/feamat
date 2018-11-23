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
elseif (n_gauss == 6)
    % coefficients taken from 
    % https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    gp = [0.816847572980459  0.091576213509771;
          0.091576213509771  0.816847572980459;
          0.091576213509771  0.091576213509771;
          0.108103018168070  0.445948490915965;
          0.445948490915965  0.108103018168070;
          0.445948490915965  0.445948490915965]';
    weights = [0.109951743655322
               0.109951743655322
               0.109951743655322
               0.223381589678011
               0.223381589678011
               0.223381589678011];
    order = 4;
elseif (n_gauss == 7)
    gp = [0.33333333333333333  0.33333333333333333;
          0.79742698535308720  0.10128650732345633;
          0.10128650732345633  0.79742698535308720;
          0.10128650732345633  0.10128650732345633;
          0.05971587178976981  0.47014206410511505;
          0.47014206410511505  0.05971587178976981;
          0.47014206410511505  0.47014206410511505]';
    weights = [0.22500000000000000
               0.12593918054482717
               0.12593918054482717
               0.12593918054482717
               0.13239415278850616
               0.13239415278850616
               0.13239415278850616];
    order = 5;
else
   error(['No 2D Gauss integration rule with ',num2str(n_gauss)',' points!!']); 
end


end

