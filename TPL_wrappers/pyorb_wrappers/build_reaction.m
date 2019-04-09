function c = build_reaction( param, current_model )

    if (strcmp( current_model, 'thermal_block' ) || strcmp( current_model, 'nonaffine_thermal_block') || strcmp( current_model, 'nonaffine' ) )
        c = @(x) 0*x(1,:) + 0*x(2,:);
    end
   
end