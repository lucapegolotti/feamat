function [sol] = solve_parameter( param, fem_specifics )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
% output=
%           sol: struct containing the solution

    [~, fespace] = set_fem_simulation( fem_specifics );

    dirichlet_functions = @(x) [0;0;0;0];
    neumann_functions = @(x) [1;0;0;0];

    % forcing term
    f = @(x) 0*x(1,:);55

    current_model = fem_specifics.model;

    if strcmp( current_model, 'nonaffine_thermal_block' )
       f = @(x) ( 1. / param(6) ) ... 
           * exp( - ( ( x(1,:)-param(4) ) .* ( x(1,:)-param(4) ) + ( x(2,:)-param(5) ) .* ( x(2,:)-param(5) ) ) / param(6) );
    end
    
    if strcmp( current_model, 'thermal_block' ) || strcmp( current_model, 'nonaffine_thermal_block' )
        mu = @(x) param(1)*(x(1,:)<0.5).*(x(2,:)<0.5) ...
        + param(2)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
        + param(3)*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
        + 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
    end
    
    if strcmp( current_model, 'nonaffine' )
        mu = @(x) ( 1. / param(3) ) ... 
           * exp( - ( ( x(1,:)-param(1) ) .* ( x(1,:)-param(1) ) + ( x(2,:)-param(2) ) .* ( x(2,:)-param(2) ) ) / param(3) );
    end

    [ A, b ] = assembler_poisson( fespace,f,mu,dirichlet_functions,neumann_functions );

    u  = A \ b;

    sol.u = u;

end

