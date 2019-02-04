function [array] = assemble_fom_rhs( param, fem_specifics, varargin )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
% output=
%           array: struct containing the stiffness matrix in COO format

    [~, fespace] = set_fem_simulation( fem_specifics, fem_specifics.mesh_name );
    
    % forcing term
    f = @(x) 0*x(1,:);

    current_model = fem_specifics.model;
    current_dirichlet = fem_specifics.use_nonhomogeneous_dirichlet;

    mu = build_diffusion( param, current_model );
    
    if strcmp( current_model, 'nonaffine' )

        f = @(x) 0*x(1,:) + 1;
        
        dirichlet_functions = @(x) [0;0;0;0];
        neumann_functions   = @(x) [0;0;0;0];
        
    end

    if strcmp( current_dirichlet, 'Y' )
        non_hom_dirichlet_functions = @(x) [1;0;0;0];
        if nargin == 2
            [ A, b ] = assembler_poisson( fespace, f, mu, dirichlet_functions, neumann_functions );
            uL = b * 0.0;
            uL = apply_dirichlet_bc_rhs( uL, fespace, non_hom_dirichlet_functions );
            b = b - A * uL;
            b = apply_dirichlet_bc_rhs( b, fespace, dirichlet_functions );
        else
%           varargin{1} is the element list;
%           varargin{2} are the indices;
            uL = zeros( size(fespace.nodes, 1), 1 );
            uL = apply_dirichlet_bc_rhs( uL, fespace, non_hom_dirichlet_functions );
            b = assemble_rhs_elementlist( fespace, f, mu, varargin{1}, uL );
            b = apply_dirichlet_bc_rhs( b, fespace, dirichlet_functions );
        end
    else
        [ ~, b ] = assembler_poisson( fespace, f, mu, dirichlet_functions, neumann_functions );
    end

    if nargin > 2
        b = b(varargin{2});
    end

    array.f = b;
end

