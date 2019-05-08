function [array] = build_fom_affine_components( operator, fem_specifics )
% Assemble fom affine matrix for elliptic scalar problems
% input=
%           operator: operator corresponding to stiffness matrix or RHS 
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           array: struct containing the affine stiffness matrices in COO format
    
    considered_model = fem_specifics.model;
    
    if ( strcmp( considered_model, 'thermal_block' ) == 0 ) && ( strcmp( operator, 'f' ) == 0 )
        error('This operator for the chosen model is not supported');
    end
    
    if isfield(fem_specifics,'final_time')
        [~, fespace] = set_fem_simulation( fem_specifics,[1;1;0;1] );
    else
        [~, fespace] = set_fem_simulation( fem_specifics );
    end

    dirichlet_functions = @(x) [0;0;0;0];
    if (isfield(fem_specifics,'final_time'))
        neumann_functions = @(x) [0;0;0;0];
    else
        neumann_functions = @(x) [1;0;0;0];
    end

    % forcing term
    if (isfield(fem_specifics,'final_time'))
        f = @(x,t) 1 * exp(-t.^2) + 0.*x(1,:);
    else
        f = @(x) 0*x(1,:); 
    end

    current_model = fem_specifics.model;

    if strcmp( current_model, 'nonaffine' )
        f = @(x) 0*x(1,:) + 1;
        dirichlet_functions = @(x) [0;0;0;0];
        neumann_functions   = @(x) [0;0;0;0];
    end
    
    if operator == 'A'
        
        if (isfield(fem_specifics,'final_time'))
            temp_f = @(x) f(x,0);
        else
            temp_f = f;
        end

        mu = @(x) (x(1,:)<0.5).*(x(2,:)<0.5);
        [ A, ~ ] = assembler_poisson( fespace,temp_f,mu,dirichlet_functions,neumann_functions );
        [i,j,val] = find( A );
        array.A0 = [i,j,val];        

        mu = @(x) (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5);
        [ A, ~ ] = assembler_poisson( fespace,temp_f,mu,dirichlet_functions,neumann_functions );
        [i,j,val] = find( A );
        array.A1 = [i,j,val];       

        mu = @(x) (x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0);
        [ A, ~ ] = assembler_poisson( fespace,temp_f,mu,dirichlet_functions,neumann_functions );
        [i,j,val] = find( A );
        array.A2 = [i,j,val];       

        mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
        [ A, ~ ] = assembler_poisson( fespace,temp_f,mu,dirichlet_functions,neumann_functions );
        [i,j,val] = find( A );
        array.A3 = [i,j,val];
    end
    
    if (isfield(fem_specifics, 'final_time'))
        if operator == 'M'
            M = assemble_mass(fespace);
            [i,j,val] = find( M );
            array.M0 = [i,j,val];
        end
    end
    
    if operator == 'f'
        
        if(isfield(fem_specifics,'final_time'))
            
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
 
            mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0);
            [ ~, b ] = assembler_poisson( fespace,f,mu,dirichlet_functions,neumann_functions );
                
        end
        
        array.f0 = b;
            
     end
        

end
    



