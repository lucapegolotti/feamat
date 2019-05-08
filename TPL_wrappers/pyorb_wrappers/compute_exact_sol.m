function [tt,exact_sol] = compute_exact_sol( param, fem_specifics, bc_flags, dirichlet_functions, neumann_functions, f_s, f_t, u_init )
% Assemble fom matrix for elliptic unsteady scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the time marching scheme
% output=
%           exact_sol: exact solution computed using ode23t

    [~, fespace] = set_fem_simulation( fem_specifics, bc_flags );
    
    T = fem_specifics.final_time;
    
    dt = T / 10000;

    current_model = fem_specifics.model;

%     if strcmp( current_model, 'nonaffine_thermal_block' )
%        f_s = @(x) ( 1. / param(6) ) ... 
%            * exp( - ( ( x(1,:)-param(4) ) .* ( x(1,:)-param(4) ) + ( x(2,:)-param(5) ) .* ( x(2,:)-param(5) ) ) / param(6) );
%     end
    
    mu_x = build_diffusion( param, current_model );
    
    c_x = build_reaction( param, current_model );
    
%     if strcmp( current_model, 'nonaffine' )        
%         f_s = @(x) 0*x(1,:) + 1;
%         dirichlet_functions = @(x) [0;0;0;0];
%         neumann_functions   = @(x) [0;0;0;0];
%     end

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
    u0 = u_init(fespace.nodes(:,1:2));  
    
    % definition of the time interval 
    tspan = 0:dt:fem_specifics.final_time;
    
    % include mass matrix option
    opts = odeset('Mass',M_time, 'RelTol', 1e-8, 'AbsTol',1e-8);
    
    % adjust rhs vector to take care of the time dependent part of the
    % forcing term
    b_adj = adjust_rhs(b,f_t,fespace);
    
    % resolution
    [tt, exact_sol] = ode23t(@(t,u) -(A+M)*u + b_adj(t).*f_t(t),tspan, u0, opts);
end   
    
    