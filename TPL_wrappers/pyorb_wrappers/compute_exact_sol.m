function [reference_sol] = compute_exact_sol( param, fem_specifics, bc_flags, dirichlet_functions, ...
                                                                        neumann_functions, f_s, f_t, u_init, timestep_number )
% Computing the "exact" solution of the problem using ode23t
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the time marching scheme
%           bc_flags: flags encoding the boundary conditions of the problem
%           dirichlet_functions: dirichlet boundary conditions
%           neumann_functions: neumann boundary conditions
%           f_s: spatial component of the forcing term
%           f_t: time component of the forcing term
%           u_init: initial condition
%           timestep_number: number of timesteps at which the exact
%           solution has to be computed
% output=
%           sol: struct containing the computed solution at the desired
%           timesteps and the timesteps themselves


    [~, fespace] = set_fem_simulation( fem_specifics, bc_flags );
    
    if timestep_number ~= -999.0
        T =  cast((fem_specifics.final_time / ...
                        fem_specifics.number_of_time_instances) * timestep_number, 'double');
         % divinding the reference timestep into 10-times smaller timesteps to
         % gain accuracy
         %dt = T / (10*timestep_number);
         dt = T / timestep_number;
    else
        T = cast(fem_specifics.final_time, 'double');
        dt = T / 1000;
    end

    current_model = fem_specifics.model;
    
    mu_x = build_diffusion( param, current_model );
    
    c_x = build_reaction( param, current_model );

    % evaluation of boundary conditions
    bc_flags = fespace.bc;
    thereisneumann = 1;
    if (length(find(bc_flags)) == length(bc_flags))
        thereisneumann = 0;
    end
    
    % evaluation of the time independent matrices and rhs
    A_no_bc = assemble_stiffness( mu_x, fespace );
    M_no_bc = assemble_mass( fespace, c_x );
    M_no_bc_time = assemble_mass( fespace );
    b_no_bc = assemble_rhs( fespace, f_s );
    
    % assembling of the lhs matrix, considering bc
    A = apply_dirichlet_bc_matrix(A_no_bc,fespace,1);
    M = apply_dirichlet_bc_matrix(M_no_bc,fespace,1);
    M_time = apply_dirichlet_bc_matrix(M_no_bc_time,fespace,0);
    
    % applying time independent bc on rhs
    if (thereisneumann)
        b = apply_neumann_bc(b_no_bc,fespace,neumann_functions);
        b = apply_dirichlet_bc_rhs(b,fespace,dirichlet_functions);
    else 
        b = apply_dirichlet_bc_rhs(b_no_bc,fespace,dirichlet_functions);
    end
    
    % evaluation of the initial condition
    y0 = u_init(fespace.nodes(:,1:2));  
    
    % definition of the time interval 
    tspan = 0:dt:T;
    
    % include mass matrix option
    opts = odeset('Mass', M_time, 'RelTol', 1e-8, 'AbsTol',1e-8);
    
    % adjust rhs vector to take care of the time dependent part of the
    % forcing term
    %b_adj = adjust_rhs(b, f_t, fespace);
    
    % resolution
    [tt, exact_sol] = ode23t(@(t,u) -(A+M)*u + b*f_t(t), tspan, y0, opts);
    exact_sol = exact_sol';
    
     if length(tspan) == 2
        tt = tt([1,end]);
        exact_sol = exact_sol(:, [1,end]);
     end
    
    % extraction of the timesteps of interest
    if timestep_number ~= -999.0
        %reference_sol.t_exact = tt(11:10:end);
        %reference_sol.u_exact = exact_sol(:, 11:10:end);
        reference_sol.t_exact = tt(2:end);
        reference_sol.u_exact = exact_sol(:, 2:end);
    else
        reference_sol.t_exact = tt;
        reference_sol.u_exact = exact_sol;
    end
    
end   
    
    