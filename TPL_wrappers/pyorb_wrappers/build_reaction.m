function c = build_reaction( param, current_model )
% Assemble the reaction coefficient as a function of space
% input=
%           param: value of the characteristic parameters
%           current_model: model for which the dreaction has to be
%           assembled as a function of space. It can be either
%           'thermal_block' or 'non_affine_thermal_block' 
% output=
%           c: space-dependent reaction coefficient

    if (strcmp( current_model, 'thermal_block' ) || strcmp( current_model, 'nonaffine_thermal_block') || strcmp( current_model, 'nonaffine' ) )
        c = @(x) 0*x(1,:) + 0*x(2,:);
    end
   
end