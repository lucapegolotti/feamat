function [sol] = get_exact_sol( params, fem_specifics, timestep_number, varargin)
% Get the "exact" solution of the problem, via ode23t. Available only for
% thermal block unstaedy problem, over the specific test cases that have
% been considered in pyorb
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the time marching scheme
%           timestep_number: number of initial timesteps where the exact
%           solution has to be evaluated
%           varargin: test case number (optional)
% output=
%           sol: struct containing the exact solution at the desired
%           timesteps and the timesteps themselves

%dbstop  in compute_exact_sol at 67

if isfield(fem_specifics, 'final_time') && strcmp(fem_specifics.model, 'thermal_block')

    if nargin > 2
            caso = varargin{1};
        else
            caso = 1;
    end
    
    timestep_number = cast(timestep_number, 'double');

    switch caso
        case 1
            
            bc_flags = [1;1;0;1];
            dirichlet_functions = @(x)   [0;0;0;0];
            neumann_functions   = @(x)   [0;0;0;0];

            % space forcing term
            f_s = @(x) 1 + 0*x(1,:) + 0*x(2,:);

            % time forcing term
            f_t = @(t) exp(-t.^2);

            % initial condition
            u_init = @(x) 0*x(:,1) + 0*x(:,2);
            
        case 2
            
            bc_flags = [1;1;1;1];
            dirichlet_functions = @(x)   [0;0;0;0];
            neumann_functions   = @(x)   [0;0;0;0];

            % space forcing term
            f_s = @(x) -(x(1,:)-x(1,:).^2).*(x(2,:)-x(2,:).^2) + 2*(x(1,:)+x(2,:)-x(1,:).^2-x(2,:).^2);

            % time forcing term
            f_t = @(t) exp(-t);

            % initial condition
            u_init = @(x) (x(:,1)-x(:,1).^2).*(x(:,2)-x(:,2).^2);
            
    end

   %computation of the exact solution
   sol = compute_exact_sol(params, fem_specifics, bc_flags, dirichlet_functions, neumann_functions, f_s, f_t, u_init, timestep_number);
    
else
    
    %exact_solution is not defined if the problem is different from an unsteady thermal block
    sol.t_exact = nan;
    sol.u_exact = nan; 

end

           
end