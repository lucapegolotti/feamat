function [gp,weights,n_gauss] = gauss_points( order )

if (order == 2)
    gp = [1/6 1/6; 2/3 1/6; 1/6 2/3]';
    weights = [1/3 1/3 1/3];
    n_gauss = 3;
else
   error(['Order ', num2str(order),' of gauss integration not implemented!']); 
end


end

