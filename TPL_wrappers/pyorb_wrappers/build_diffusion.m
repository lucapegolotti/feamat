function mu = build_diffusion( param, current_model )
% Assemble the diffusion coefficient as a function of space
% input=
%           param: value of the characteristic parameters
%           current_model: model for which the diffusion has to be
%           assembled as a function of space. It can be either
%           'thermal_block' or 'non_affine_thermal_block' or 'non_affine'
% output=
%           mu: space-dependent diffusion coefficient

    if strcmp( current_model, 'thermal_block' ) || strcmp( current_model, 'nonaffine_thermal_block' )
        mu = @(x) param(1)*(x(1,:)<0.5).*(x(2,:)<0.5) ...
        + param(2)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
        + param(3)*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
        + 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
    end
    
    if strcmp( current_model, 'nonaffine' )
        mu = @(x) 1. * ( param(3) + ( 1. / param(3) ) ... 
           * exp( - ( ( x(1,:)-param(1) ) .* ( x(1,:)-param(1) ) ...
                    + ( x(2,:)-param(2) ) .* ( x(2,:)-param(2) ) ) / param(3) ) );
                
    end

end