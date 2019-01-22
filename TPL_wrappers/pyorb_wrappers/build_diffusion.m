function mu = build_diffusion( param, current_model )

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