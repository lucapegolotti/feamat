function [gp,weights,n_gauss] = gauss_points1D( order )

if (order == 3)
    gp = [-sqrt(1/3); sqrt(1/3)]';
    weights = [1 1];
    n_gauss = 2;
elseif (order == 5)
    gp = [-sqrt(3/5);0;sqrt(3/5)]';
    weights = [5/9 8/9 5/9];
    n_gauss = 3;
else
   error(['Order ', num2str(order),' of gauss integration 1D not implemented!']); 
end


end



