clearvars
close all
clc 

%% Convergence test

fem_specifics.number_of_elements = 100;
fem_specifics.polynomial_degree = 'P1';
fem_specifics.model = 'thermal_block';
fem_specifics.use_nonhomogeneous_dirichlet = 'N';
fem_specifics.mesh_name = 'cube100x100';

% New fields in time dependent case!!!!
fem_specifics.final_time = 0.1;
fem_specifics.theta = 1.0;

params = [1.00, 1.00, 1.00]; 

bc_flags = [1;1;0;1];
dirichlet_functions = @(x)   [0;0;0;0];
neumann_functions   = @(x)   [0;0;0;0];

% space forcing term
%f_s = @(x) (8*pi^2-1) * sin(2*pi*x(1,:)) .* sin(2*pi*x(2,:));
f_s = @(x) 1 + 0*x(1,:) + 0*x(2,:);
    
% time forcing term
f_t = @(t) exp(-t.^2);
    
% initial condition
%u_init = @(x) sin(2*pi*x(:,1)) .* sin(2*pi*x(:,2));
u_init = @(x) 0*x(:,1) + 0*x(:,2);
    

%times = [ 1 2 4 8 16 32 64 128 256 ]';
times = [ 5 10 20 40 80 160 320 ]';
time_steps = fem_specifics.final_time ./ times;

[~, fespace] = set_fem_simulation( fem_specifics, bc_flags );

M_no_bc_time = assemble_mass( fespace );
M = apply_dirichlet_bc_matrix(M_no_bc_time,fespace,0);
% M = M_no_bc_time;

%%
% computation of the exact solution
% T = fem_specifics.final_time;
% exact_sol = @(x) exp(-T) * sin(2*pi*x(1)) .* sin(2*pi*x(2));
% grad_exact_sol = @(x) [2*pi*cos(2*pi*x(1))*sin(2*pi*x(2)); 2*pi*sin(2*pi*x(1))*cos(2*pi*x(2))] * exp(-T);  

[t_exact, exact_sol] = compute_exact_sol(params, fem_specifics, bc_flags, dirichlet_functions, neumann_functions, f_s, f_t, u_init);


%%
% computation of numerical solutions for different time steps and errors
err_L2 = zeros(length(times),1);
err_H1 = zeros(length(times),1);

err_L2_new = zeros(length(times),1);

for count = 1:length(times)
    
    fem_specifics.number_of_time_instances = times(count);
    sol = solve_parameter( params, fem_specifics );
    %err_L2(count) = norm(sol.u(:,end) - exact_sol(end,:)',2);
    err_L2(count) = compute_norm(fespace,sol.u(:,end) - exact_sol(end,:)','L2');
    err_H1(count) = compute_norm(fespace,sol.u(:,end) - exact_sol(end,:)','H1');

    u_end = sol.u(:,end);
    u_ex_end = exact_sol(end,:)';
    err_end = u_end - u_ex_end;
    err_L2_new(count) = sqrt(err_end'*M*err_end);
    
    %err_L2(count) = compute_L2_error(fespace,sol.u(:,end),exact_sol(:,end));
    %err_H1(count) = compute_H1_error(fespace,sol.u(:,end),exact_sol, grad_exact_sol);
    
end

figure
loglog(time_steps, err_L2, '-o');
hold on
grid on
loglog(time_steps, time_steps, '--');
title('L2 error');
legend('L2 error', 'dt');
hold off

order_L2 = log(err_L2(2:end)./err_L2(1:end-1)) ./ log(time_steps(2:end)./time_steps(1:end-1));

figure
loglog(time_steps, err_H1, '-o');
hold on
grid on
loglog(time_steps, time_steps, '--');
title('H1 error');
legend('H1 error', 'dt');

order_H1 = log(err_H1(2:end)./err_H1(1:end-1)) ./ log(time_steps(2:end)./time_steps(1:end-1));

figure
loglog(time_steps, err_L2_new, '-o');
hold on
grid on
loglog(time_steps, time_steps, '--');
title('L2 error new');
legend('L2 error new', 'dt');

order_L2_new = log(err_L2_new(2:end)./err_L2_new(1:end-1)) ./ log(time_steps(2:end)./time_steps(1:end-1));

%% Plot Exact solution
x = fespace.nodes(:,1);
y = fespace.nodes(:,2);
%exact_sol = exp(-T) * sin(2*pi*x) .* sin(2*pi*y);

figure
for i=1:length(t_exact)
    plot3(x,y,exact_sol(i,:));
    grid on
    title('Reference solution');
    zlim([0;0.2]);
    pause(0.001)
end

%% Plot numerical solution

figure
for i=1:fem_specifics.number_of_time_instances
    plot3(x,y,sol.u(:,i));
    grid on
    title('Numerical solution');
    zlim([0;0.2]);
    pause(0.001)
end

