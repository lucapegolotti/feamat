function [array] = assemble_fom_rhs( param, fem_specifics, varargin )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
% output=
%           array: struct containing the stiffness matrix in COO format

    [~, fespace] = set_fem_simulation( fem_specifics );

    % forcing term
    if (isfield(fem_specifics,'final_time'))
        f = @(x,t) 1 * exp(-t.^2) + 0.*x(1,:);
    else
        f = @(x) 0*x(1,:);
    end

    current_model = fem_specifics.model;
    current_dirichlet = fem_specifics.use_nonhomogeneous_dirichlet;

    mu = build_diffusion( param, current_model );
    
    if strcmp( current_model, 'nonaffine' )

        f = @(x) 0*x(1,:) + 1;
        
        dirichlet_functions = @(x) [0;0;0;0];
        neumann_functions   = @(x) [0;0;0;0];
        
    end
    
    if (isfield(fem_specifics,'final_time'))
        dim = (fem_specifics.number_of_elements+1)^2;
        b = zeros(dim * fem_specifics.number_of_time_instances,1);
        dt = fem_specifics.final_time / fem_specifics.number_of_time_instances;

        index = 0;
        for time = dt:dt:fem_specifics.final_time

             temp_f = @(x) f(x,time);

             mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0);
             [ ~, b(index*dim+1:(index+1)*dim) ] = assembler_poisson( fespace,temp_f,mu,dirichlet_functions,neumann_functions );

             index = index + 1;

        end
        
    else
        
        if strcmp( current_dirichlet, 'Y' )
            non_hom_dirichlet_functions = @(x) [1;0;0;0];
            [ A, b ] = assembler_poisson( fespace, f, mu, non_hom_dirichlet_functions, neumann_functions );
            uL = b * 0.0;
            uL = apply_dirichlet_bc_rhs( uL, fespace, non_hom_dirichlet_functions );
            b = b - A * uL;
            b = apply_dirichlet_bc_rhs( b, fespace, dirichlet_functions );
        else
        [ ~, b ] = assembler_poisson( fespace, f, mu, dirichlet_functions, neumann_functions );
        end

        if nargin > 2
            b = b(varargin{2});
        end

    end
    
    array.f = b;
    
end

