function [A] = assemble_diffusive_term(fespace, mu)
    % Assemble the stiffness matrix associated to a scalar/tensorial,
    % uniform/space-dependent diffusivity coefficient mu.
    % Input:
    %           fespace: finite elements space object
    %           mu: diffusivity coefficient. It is a space-dependent or 
    %           constant 2x2 matrix.
    % Output:   diffusivity matrix.

    if(~isa(mu,'function_handle'))
        A = assemble_diffusive_term_component(fespace, 1, 1, mu(1,1)) + ...
            assemble_diffusive_term_component(fespace, 2, 2, mu(2,2)) + ...
            assemble_diffusive_term_component(fespace, 1, 2, mu(1,2)) + ...
            assemble_diffusive_term_component(fespace, 2, 1, mu(2,1));
    else
        e1 = [1;0];
        e2 = [0;1];
        mu11 = @(x) e1'*(mu(x)*e1);
        mu22 = @(x) e2'*(mu(x)*e2);
        mu12 = @(x) e1'*(mu(x)*e2);
        mu21 = @(x) e2'*(mu(x)*e1);
        A = assemble_diffusive_term_component(fespace, 1, 1, mu11) + ...
            assemble_diffusive_term_component(fespace, 2, 2, mu22) + ...
            assemble_diffusive_term_component(fespace, 1, 2, mu12) + ...
            assemble_diffusive_term_component(fespace, 2, 1, mu21);
    end
end