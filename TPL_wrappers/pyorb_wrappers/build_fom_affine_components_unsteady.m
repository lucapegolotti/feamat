function [array] = build_fom_affine_components_unsteady( operator, fem_specifics, varargin )
% Assemble fom affine matrix for elliptic scalar problems
% input=
%           operator: operator corresponding to stiffness matrix or RHS 
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model         
%           varargin: if passed, it contains the reference test number as
%           second argument (used in thermal block unsteady!) and the
%           desired number of time instances as first argument
% output=
%           array: struct containing the affine stiffness matrices in COO format

    
    considered_model = fem_specifics.model;
    
    if isfield(fem_specifics, 'final_time') && isfield(fem_specifics, 'number_of_time_instances')
        goOn = 1;
    else
        disp("The problem does not appear to be unsteady. The steady version of this function is called")
        array = build_fom_affine_components(operator, fem_specifics);
        goOn = 0;
    end
    
    if goOn > 0
        if ( strcmp( considered_model, 'thermal_block' ) == 0 ) && ( strcmp( operator, 'f' ) == 0 )
            error('This operator for the chosen model is not supported');
        end

        [~, fespace] = set_fem_simulation( fem_specifics, [1;1;0;1] );

        dirichlet_functions = @(x) [0;0;0;0];
        neumann_functions = @(x) [0;0;0;0];

        % forcing term
        if  (strcmp(considered_model, 'thermal_block')==1)
            if nargin <= 4 || (nargin == 4 && varargin{2}==1)
                f_space = @(x) 1 + 0.*x(1,:) + 0.*x(2,:);
                f_time = @(t) exp(-t.^2) ;
                f = @(x,t) f_space(x) .* f_time(t);
            elseif nargin==4 && varargin{2}==2
                f_space = @(x) 1 + 0.*x(1,:) + 0.*x(2,:);
                f_time = @(t) exp(-t) ;
                f = @(x,t) f_space(x) .* f_time(t);
            else
                f_space = @(x) 1 + 0.*x(1,:) + 0.*x(2,:);
                f_time = @(t) 0.*t;
                f = @(x,t) f_space(x) .* f_time(t);
            end
        else
             f_space = @(x) 1 + 0.*x(1,:) + 0.*x(2,:);
             f_time = @(t) 0.*t;
             f = @(x,t) f_space(x) .* f_time(t);
        end

        if strcmp(operator, 'A')

            mu = @(x) (x(1,:)<0.5).*(x(2,:)<0.5);
            [ A, ~ ] = assembler_poisson( fespace,f_space,mu,dirichlet_functions,neumann_functions );
            [i,j,val] = find( A );
            array.A0 = [i,j,val];        

            mu = @(x) (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5);
            [ A, ~ ] = assembler_poisson( fespace,f_space,mu,dirichlet_functions,neumann_functions );
            [i,j,val] = find( A );
            array.A1 = [i,j,val];       

            mu = @(x) (x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0);
            [ A, ~ ] = assembler_poisson( fespace,f_space,mu,dirichlet_functions,neumann_functions );
            [i,j,val] = find( A );
            array.A2 = [i,j,val];       

            mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
            [ A, ~ ] = assembler_poisson( fespace,f_space,mu,dirichlet_functions,neumann_functions );
            [i,j,val] = find( A );
            array.A3 = [i,j,val];
        end
        

        if strcmp(operator, 'M')
            M = assemble_mass(fespace);
            [i,j,val] = find( M );
            array.M0 = [i,j,val];
        end
        

        if strcmp(operator, 'f')

            dim = (fem_specifics.number_of_elements+1)^2;
               
            if nargin == 2
                b = zeros(dim * fem_specifics.number_of_time_instances, 1);
                dt = fem_specifics.final_time / fem_specifics.number_of_time_instances;
            else
                b = zeros(dim * varargin{1}, 1);
                dt = fem_specifics.final_time / varargin{1};
            end

            index = 0;
            for time = 0:dt:fem_specifics.final_time

                temp_f = @(x) f(x,time);

                mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0);
                [ ~, b(index*dim+1:(index+1)*dim) ] = assembler_poisson(...
                                                        fespace,temp_f,mu,dirichlet_functions,neumann_functions );

                index = index + 1;

            end

            array.f0 = b;

        end
        
        if strcmp(operator, 'f_space')
            
            mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
             [ ~, b ] = assembler_poisson( fespace,f_space,mu,dirichlet_functions,neumann_functions );
             array.f_space0 = b;
             
        end
                 
    end
        
end
    



