function [sol] = solve_parameter_unsteady( param, fem_specifics, varargin )
%Solve the unsteady fom problem using a suitable multistep time-marching
%scheme
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the time marching scheme
%           varargin: test case number (optional)
% output=
%           sol: struct containing the solution

    if nargin > 2
        caso = varargin{1};
    else
        caso = 1;
    end
    
    switch caso
        case 1
            bc = [1;1;0;1];
        case 2
            bc = [1;1;1;1];
    end
    [~, fespace] = set_fem_simulation( fem_specifics, bc );
    
    % TERMS TO BE ADDED IN fem_specifics!!!!
    n_time_instances = cast(fem_specifics.number_of_time_instances, 'double');
    T = cast(fem_specifics.final_time, 'double');
    Theta = fem_specifics.theta; % parameter of theta-method in [0,1]
    step_number = cast(fem_specifics.step_number_fom, 'double'); %number of steps in multistep implementation
    
    % computation of the time step
    dt = cast(T / n_time_instances, 'double');
    
     switch caso
        case 1

            dirichlet_functions = @(x)   [0;0;0;0];
            neumann_functions   = @(x)   [0;0;0;0];

            % space forcing term
            f_s = @(x) 1 + 0*x(1,:) + 0*x(2,:);

            % time forcing term
            f_t = @(t) exp(-t.^2);

            % initial condition
            u_init = @(x) 0*x(:,1) + 0*x(:,2);
            
        case 2
            
            dirichlet_functions = @(x)   [0;0;0;0];
            neumann_functions   = @(x)   [0;0;0;0];

            % space forcing term
            f_s = @(x) -(x(1,:)-x(1,:).^2).*(x(2,:)-x(2,:).^2) + 2*(x(1,:)+x(2,:)-x(1,:).^2-x(2,:).^2);

            % time forcing term
            f_t = @(t) exp(-t);

            % initial condition
            u_init = @(x) (x(:,1)-x(:,1).^2).*(x(:,2)-x(:,2).^2);
    end

    current_model = fem_specifics.model;
    
    mu_x = build_diffusion( param, current_model );
    
    c_x = build_reaction( param, current_model );
    
    % evaluation of the time independent matrices and rhs
    A_no_bc = assemble_stiffness( mu_x, fespace );
    M_no_bc = assemble_mass( fespace, c_x );
    M_no_bc_time = assemble_mass( fespace );
    b_no_bc = assemble_rhs( fespace, f_s );
    
    % getting the lhs multiplicative coefficient
    if step_number == 1
        coeffs = get_multistep_coefficients(step_number, Theta);
    else
        coeffs = get_multistep_coefficients(step_number);
    end
    
    % assembling of the lhs matrix, considering bc
    LHS = M_no_bc_time + coeffs(1) * dt * (A_no_bc + M_no_bc );
    LHS = apply_dirichlet_bc_matrix(LHS, fespace, 1);
    
    % performing LU fatorization of the matrix outside of the time loop
    [L,U,P] = lu(LHS);
    
    % evaluation of boundary conditions
    bc_flags = fespace.bc;
    thereisneumann = 1;
    if (length(find(bc_flags)) == length(bc_flags))
        thereisneumann = 0;
    end
    
    % defining the fem dimension
    fem_dimension = length(fespace.nodes(:,1));
    
    % initialization of the solution matrix
    u = zeros(fem_dimension, n_time_instances);

    % initial time
    t = 0.0;
    
    % evaluation of the initial condition
    u(:,1) = u_init(fespace.nodes(:,1:2));
    
    %evaluation of the initial timesteps, if needed by the multistep method
    if step_number >= 2
         sol = compute_exact_sol(param, fem_specifics, bc_flags, dirichlet_functions, ...
                  neumann_functions, f_s, f_t, u_init, step_number-1);
         u(:, 2:step_number) = sol.u_exact;       
    end
     
     %pre-assembling of the vectors that will be involved in the rhs
     %computation
     if (step_number > 1) || ((step_number == 1) && (coeffs(end) > 0))
         rhs_mat = zeros(fem_dimension, step_number);
         for step = 1:step_number
             rhs_mat(:, step) = (A_no_bc + M_no_bc) * u(:, step);
         end
     end
 
    % evaluation of the time part of the force at the needed initial time
    % instants 
    f0 = zeros(step_number, 1);
    for step = 1:step_number
        f0(step) = f_t(t+(step_number-1)*dt);
    end
         
    % starting time loop
    count = step_number;
    for t = step_number*dt:dt:T   
        
        % evaluation of the time part of the force at new time instant
        f1 = f_t(t);
        
        % multiplication of the time forces with the rhs, without considering
        % bc imposition  
        F = zeros(length(fespace.nodes(:,1)), step_number+1);
        for step = 1:step_number
            F(:,step) = f0(step) * b_no_bc;
        end
        F(:,end) = f1 * b_no_bc;

        % assembling of the rhs vector, without considering bc
       RHS =  M_no_bc_time * u(:, count);
       for step = 1:step_number+1
           RHS = RHS + coeffs(step) * dt * F(:, end-step+1);
            if step > 1 && coeffs(step) > 0
                RHS = RHS - coeffs(step) * dt * rhs_mat(:, end-step+2);
            end
       end      
        
        % imposition of bc in classical way rhs
        if (thereisneumann)
           RHS = apply_neumann_bc(RHS, fespace, neumann_functions); 
        end

        RHS = apply_dirichlet_bc_rhs(RHS, fespace, dirichlet_functions);

        % system resolution using LU factorization
        temp = L \ (P*RHS);
        u(:,count+1) = U \ temp; 

        u(:,count+1) = apply_dirichlet_bc_rhs(u(:,count+1), fespace, dirichlet_functions);

        % update the time part of the force and the counter  
        if step_number >= 2
            f0(1:end-1) = f0(2:end);
            rhs_mat(:, 1:end-1) = rhs_mat(:, 2:end);
        end
        f0(end) = f1;
        if (step_number>1) || (step_number==1 && coeffs(end)>0)
            rhs_mat(:, end) = (A_no_bc + M_no_bc) * u(:, count+1);
        end
        count = count + 1;
        
    end
 
    sol.u = u;

end

