function [sol] = solve_parameter( param, fem_specifics, varargin )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
%           varargin: test case number (optional)
% output=
%           sol: struct containing the solution

    tic;

    if (isfield( fem_specifics, 'final_time' ))
        % in this case it means that the problem is unsteady, since more
        % fields (like 'final_time') have been added in such a case!!!
        if nargin > 2
            sol = solve_parameter_unsteady( param, fem_specifics, varargin{1} );
        else
            sol = solve_parameter_unsteady( param, fem_specifics );
        end
    else
        % otherwise we proceed in the usual fashion, with the steady
        % problem
        
        [~, fespace] = set_fem_simulation( fem_specifics );

        dirichlet_functions = @(x) [0;0;0;0];
        neumann_functions = @(x) [1;0;0;0];

        % forcing term
        f = @(x) 0*x(1,:);

        current_model = fem_specifics.model;
        current_dirichlet = fem_specifics.use_nonhomogeneous_dirichlet;

        if strcmp( current_model, 'nonaffine_thermal_block' )
           f = @(x) ( 1. / param(6) ) ... 
               * exp( - ( ( x(1,:)-param(4) ) .* ( x(1,:)-param(4) ) + ( x(2,:)-param(5) ) .* ( x(2,:)-param(5) ) ) / param(6) );
        end

        mu = build_diffusion( param, current_model );

        if strcmp( current_model, 'nonaffine' )        
            f = @(x) 0*x(1,:) + 1;
            dirichlet_functions = @(x) [0;0;0;0];
            neumann_functions   = @(x) [0;0;0;0];
        end

        if strcmp( current_dirichlet, 'Y' )
            non_hom_dirichlet_functions = @(x) [1;0;0;0];
            [ A, b ] = assembler_poisson( fespace, f, mu, non_hom_dirichlet_functions, neumann_functions );
            uL = b * 0.0;
            uL = apply_dirichlet_bc_rhs( uL, fespace, non_hom_dirichlet_functions );
            b = b - A * uL;
            b = apply_dirichlet_bc_rhs( b, fespace, dirichlet_functions );
            [ A, ~ ] = assembler_poisson( fespace, f, mu, dirichlet_functions, neumann_functions );
        else
            [ A, b ] = assembler_poisson( fespace, f, mu, dirichlet_functions, neumann_functions );
        end

        u  = A \ b;

        sol.u = u;
        sol.execution_time = toc;
    end
end

