clearvars
close all
clc 

%% Convergence test

spatial_dim = 50; % choose the dimension of the FE space
N_dims = 1;

% choose between case 1 (reference solution computed via ode23t) and case 2
% (reference solution computed from an analytical expression)
caso = 1; 

count_dims = 1;

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

A = assemble_stiffness(1, fespace);




