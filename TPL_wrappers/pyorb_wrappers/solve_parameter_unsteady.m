function [sol] = solve_parameter_unsteady( param, fem_specifics )
% Assemble fom matrix for elliptic unsteady scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the time marching scheme
% output=
%           sol: struct containing the solution

    [~, fespace] = set_fem_simulation( fem_specifics, [1,1,1,1] );
    
    % TERMS TO BE ADDED IN fem_specifics!!!!
    n_time_instances = fem_specifics.number_of_time_instances;
    T = fem_specifics.final_time;
    Theta = fem_specifics.theta; % parameter of theta method in (0,1) 
    
    % computation of the time step
    dt = T / n_time_instances;

    dirichlet_functions = @(x) [0;0;0;0];
    neumann_functions = @(x) [0;0;0;0];

    % space forcing term
    f_s = @(x) (8*pi^2-1) * sin(2*pi*x(1,:)) .* sin(2*pi*x(2,:));
    
    % time forcing term
    f_t = @(t) exp(-t);
    
    % initial condition
    u_init = @(x) sin(2*pi*x(:,1)) .* sin(2*pi*x(:,2));

    current_model = fem_specifics.model;

%     if strcmp( current_model, 'nonaffine_thermal_block' )
%        f_s = @(x) ( 1. / param(6) ) ... 
%            * exp( - ( ( x(1,:)-param(4) ) .* ( x(1,:)-param(4) ) + ( x(2,:)-param(5) ) .* ( x(2,:)-param(5) ) ) / param(6) );
%     end
    
    mu = build_diffusion( param, current_model );
    
%     if strcmp( current_model, 'nonaffine' )        
%         f_s = @(x) 0*x(1,:) + 1;
%         dirichlet_functions = @(x) [0;0;0;0];
%         neumann_functions   = @(x) [0;0;0;0];
%     end
    
    % evaluation of the time independent matrices and rhs
    A_no_bc = assemble_stiffness( mu, fespace );
    M_no_bc = assemble_mass( fespace );
    b_no_bc = assemble_rhs( fespace, f_s );
    
    % evaluation of boundary conditions
    bc_flags = fespace.bc;
    thereisneumann = 1;
    if (length(find(bc_flags)) == length(bc_flags))
        thereisneumann = 0;
    end
    
    % initialization of the solution matrix
    u = zeros(length(fespace.nodes(:,1)), n_time_instances);

    % initial time
    t = 0.0;
    
    % evaluation of the initial condition
    u(:,1) = u_init(fespace.nodes(:,1:2));
    
    % evaluation of the time part of the force at t=0
    f0 = f_t(t);
    
    % starting time loop
    count = 1;
    for t = dt:dt:T
        % evaluation of the time part of the force at new time instant
        f1 = f_t(t);
        
        % multiplication of the time forces with the rhs, without considering
        % bc imposition
        F0 = f0 * b_no_bc;
        F1 = f1 * b_no_bc;
        
        % assembling of the lhs matrix, without considering bc
        LHS = M_no_bc + Theta * dt * A_no_bc;
        
        % assembling of the rhs vector, without considering bc
        RHS = (M_no_bc - (1 - Theta) * dt  * A_no_bc) * u(:,count) + Theta * F1 * dt + (1 - Theta) * F0 * dt;
        
        % imposition of bc in classical way on lhs and rhs
        if (thereisneumann)
           RHS = apply_neumann_bc(RHS,fespace,neumann_functions); 
        end

        [LHS,RHS] = apply_dirichlet_bc(LHS,RHS,fespace,dirichlet_functions);

        % system resolution
        u(:,count+1) = LHS\RHS;

        % update the time part of the force and the counter   
        f0 = f1;
        count = count + 1;
        
    end
 
    sol.u = u;

end

