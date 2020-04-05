function [array] = build_fom_affine_components( operator, fem_specifics )
% Assemble the FEM approximated affine components of the stiffness matrix or 
% of the RHS forcing term for steady elliptic scalar problems
% input=
%           operator: operator corresponding to stiffness matrix ('A') or
%           to the RHS ('f')
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           array: struct containing the affine stiffness matrices (in sparse COO
%           format) or the affine RHS (in full format)
    
    considered_model = fem_specifics.model;
    
    if ( strcmp( considered_model, 'thermal_block' ) == 0 ) && ( strcmp( operator, 'f' ) == 0 )
        error('This operator for the chosen model is not supported');
    end

    [~, fespace] = set_fem_simulation( fem_specifics );

    dirichlet_functions = @(x) [0;0;0;0];
    neumann_functions = @(x) [1;0;0;0];

    % forcing term
    f = @(x) 0*x(1,:);

    current_model = fem_specifics.model;

    if strcmp( current_model, 'nonaffine' )
        f = @(x) 0*x(1,:) + 1;
        dirichlet_functions = @(x) [0;0;0;0];
        neumann_functions   = @(x) [0;0;0;0];
    end
    
    if operator == 'A'

        mu = @(x) (x(1,:)<0.5).*(x(2,:)<0.5);
        [ A, ~ ] = assembler_poisson( fespace,f,mu,dirichlet_functions,neumann_functions );

        [i,j,val] = find( A );
        array.A0 = [i,j,val];

        mu = @(x) (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5);
        [ A, ~ ] = assembler_poisson( fespace,f,mu,dirichlet_functions,neumann_functions );

        [i,j,val] = find( A );
        array.A1 = [i,j,val];

        mu = @(x) (x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0);
        [ A, ~ ] = assembler_poisson( fespace,f,mu,dirichlet_functions,neumann_functions );

        [i,j,val] = find( A );
        array.A2 = [i,j,val];

        mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
        [ A, ~ ] = assembler_poisson( fespace,f,mu,dirichlet_functions,neumann_functions );

        [i,j,val] = find( A );
        array.A3 = [i,j,val];

    end
    
    if operator == 'f'
        
        mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
        [ ~, b ] = assembler_poisson( fespace,f,mu,dirichlet_functions,neumann_functions );
        array.f0 = b;

    end
    
end