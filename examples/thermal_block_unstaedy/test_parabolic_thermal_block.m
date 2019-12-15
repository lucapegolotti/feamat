clearvars
close all
clc 

%% Convergence test

dims = [50]; % choose the dimension of the FE space
N_dims = length(dims);

% choose between case 1 (reference solution computed via ode23t) and case 2
% (reference solution computed from an analytical expression)
caso = 1; 

count_dims = 1;
for spatial_dim = dims

    fem_specifics.number_of_elements = spatial_dim;
    fem_specifics.polynomial_degree = 'P1';
    fem_specifics.model = 'thermal_block';
    fem_specifics.use_nonhomogeneous_dirichlet = 'N';                                                            
    fem_specifics.mesh_name = ['cube' num2str(spatial_dim) 'x' num2str(spatial_dim)];

    % New fields in time dependent case!!!!
    switch caso
        case 1
            fem_specifics.final_time = 0.1;
        case 2
            fem_specifics.final_time = 1.0;
    end
    fem_specifics.theta = 1.0;
    fem_specifics.step_number_fom = 4;
    fem_specifics.method = 'BDF';

    params = [9.12, 0.13, 7.6]; 
    
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

            %times = [ 1 5 10 20 40 80 160 320 ]';
            times =50;
            time_steps = fem_specifics.final_time ./ times;
            
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

            times = [ 1 2 4 8 16 32 54 128 ]';
            %times = 25;
            time_steps = fem_specifics.final_time ./ times;
            
    end
    
    [~, fespace] = set_fem_simulation( fem_specifics, bc_flags );

    T = fem_specifics.final_time;

    % computation of the exact solution
    switch caso
        case 1
            timestep_number = -999; % symbolic value to get the reference solution at all timesteps
            exact_sol = compute_exact_sol(params, fem_specifics, bc_flags, ...
                               dirichlet_functions, neumann_functions, f_s, f_t, u_init, timestep_number);
            
        case 2
            
            exact_sol = @(x) exp(-T) * (x(1)-x(1).^2).*(x(2)-x(2).^2);
            grad_exact_sol = @(x) [exp(-T) * (1-2*x(1)).*(x(2)-x(2).^2); ...
                                                 exp(-T) * (1-2*x(2)).*(x(1)-x(1).^2)];
          
    end
    
    
    % computation of numerical solutions for different time steps and errors
    err_L2 = zeros(length(times), N_dims);
    err_H1 = zeros(length(times), N_dims);
    order_L2 = zeros(length(times)-1, N_dims);
    order_H1 = zeros(length(times)-1, N_dims);

    for count = 1:length(times)

        fem_specifics.number_of_time_instances = times(count);
        sol = solve_parameter( params, fem_specifics, caso );
        
        switch caso
            case 1
                err_L2(count,count_dims) = compute_norm(fespace,sol.u(:, end) - ...
                                                            exact_sol.u_exact(:, end), 'L2');
                err_H1(count,count_dims) = compute_norm(fespace,sol.u(:, end) - ...
                                                            exact_sol.u_exact(:,end), 'H1');
                
            case 2
                
                err_L2(count,count_dims) = compute_L2_error(fespace,sol.u(:,end),exact_sol);
                err_H1(count,count_dims) = compute_H1_error(fespace,sol.u(:,end),exact_sol,grad_exact_sol);
                
        end

    end

    figure(1)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.25, 0.7, 0.75]);
    subplot(N_dims,1,count_dims)
    loglog(time_steps, err_L2(:,count_dims), '-o');
    hold on
    grid on
    loglog(time_steps, time_steps/100, '--');
    title(['L2 error on mesh: ' fem_specifics.mesh_name]);
    xlabel('Number of timesteps')
    ylabel('L2 error')
    legend('L2 error', 'dt','Location','Best');
    

    order_L2(:,count_dims) = log(err_L2(2:end,count_dims)./err_L2(1:end-1,count_dims)) ./ log(time_steps(2:end)./time_steps(1:end-1));

    figure(2)
    subplot(N_dims,1,count_dims)
    loglog(time_steps, err_H1(:,count_dims), '-o');
    hold on
    grid on
    loglog(time_steps, time_steps/100, '--');
    title(['H1 error on mesh: ' fem_specifics.mesh_name]);
    xlabel('Number of timesteps')
    ylabel('H1 error')
    legend('H1 error', 'dt','Location', 'Best');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.25, 0.7, 0.75]);

    order_H1(:,count_dims) = log(err_H1(2:end,count_dims)./err_H1(1:end-1,count_dims)) ./ log(time_steps(2:end)./time_steps(1:end-1));
    
    count_dims = count_dims +1;
end

%% Plotting Area

% Plot Exact solution
x = fespace.nodes(:,1);
y = fespace.nodes(:,2);

switch caso
    case 1
        
        figure(3)
        subplot(1,2,1)
        for i=1:10:length(exact_sol.t_exact)
            plot3(x,y,exact_sol.u_exact(:,i));
            grid on
            title('Reference solution');
            zlim([0;0.2]);
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.3, 0.75, 0.6]);
            pause(0.00001)
        end
        
    case 2
        
        figure(3)
        subplot(1,2,1)
        exact_sol = @(x,y,t) exp(-t) .* (x-x.^2).*(y-y.^2);
        t_exact = linspace(0,T,1000);
        for i=1:length(t_exact)
            exact_sol_t = @(x,y) exact_sol(x,y,t_exact(i));
            plot3(x,y,exact_sol_t(x,y));
            grid on
            title('Reference solution');
            zlim([0;0.06]);
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.3, 0.75, 0.6]);
            pause(0.00001)
        end
        
end

% Plot numerical solution
x = fespace.nodes(:,1);
y = fespace.nodes(:,2);

switch caso
    case 1

        figure(3)
        subplot(1,2,2)
        for i=1:fem_specifics.number_of_time_instances
            plot3(x,y,sol.u(:,i));
            grid on
            title('Numerical solution');
            zlim([0;0.2]);
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.3, 0.75, 0.6]);
            pause(0.00001)
        end
        
    case 2
        
        figure(3)
        subplot(1,2,2)
        for i=1:fem_specifics.number_of_time_instances
            plot3(x,y,sol.u(:,i));
            grid on
            title('Numerical solution');
            zlim([0;0.06]);
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.3, 0.75, 0.6]);
            pause(0.00001)
        end
        
end

