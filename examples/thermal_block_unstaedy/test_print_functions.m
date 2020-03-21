clearvars
close all
clc 

%% Defining the FE space

dim = 50;

fem_specifics.number_of_elements = dim;
fem_specifics.polynomial_degree = 'P1';
fem_specifics.model = 'thermal_block';
fem_specifics.use_nonhomogeneous_dirichlet = 'N';                                                            
fem_specifics.mesh_name = ['cube' num2str(dim) 'x' num2str(dim)];

% New fields in time dependent case!!!!

fem_specifics.final_time = 0.1;
fem_specifics.number_of_time_instances = 1000.0;
fem_specifics.method = 'BDF';
fem_specifics.step_number_fom = 6;
fem_specifics.theta = 1.0;

bc = [1;1;0;1];

[~, fespace] = set_fem_simulation( fem_specifics, bc );

%% Importing the RB basis files
 
basis_space = importdata('./../TESTS/RB/Unsteady/tbp_multistep/BDF/Test1/data/fem_dim_50_time_dim_1000/pod_tolerance_-5/basis.txt');
%basis_time = importdata('/home/tenderin/Desktop/TESTS/Unsteady/thermal_block_test/test_heat_equation_offline_50/basis_time_50_train.txt');


%% Plotting space basis

figure()
for i = 35:39
    subplot(2,2,i-34)
    h = plot_fe_function(basis_space(:,i),fespace);
    title(['Space Basis Function ', num2str(i)]);
end


%% Plotting Time basis

figure()
for i = 1:4
    subplot(2,2,i)
    t = linspace(0,fem_specifics.final_time, fem_specifics.number_of_time_instances);
    plot(t,basis_time(:,i),'-');
    grid on
    title(['Time Basis Function ', num2str(i)]);
end

%% Importing RB solution and error

rb_solution  = importdata('/home/tenderin/Desktop/TESTS/Unsteady/thermal_block_test/test_heat_equation_offline_50/reduced_solution50_train.txt');
fom_solution = importdata('/home/tenderin/Desktop/TESTS/Unsteady/thermal_block_test/test_heat_equation_offline_50/fom_solution50_train.txt');

rb_error = abs(rb_solution - fom_solution);

%% Printing RB solution and error

figure()

subplot(1,2,1)
plot_fe_function(rb_solution(:,1),fespace);
title('RB solution');

subplot(1,2,2)
plot_fe_function(rb_error(:,1)./fom_solution(:,1),fespace);
title('Logarithm of the error');

%% Printing video errors

figure()

for i = 1:length(rb_solution(1,:))
    subplot(1,2,1)
    plot_fe_function(rb_solution(:,i),fespace,'surf',[0 0.05]);
    title('RB solution');

    subplot(1,2,2)
    plot_fe_function(rb_error(:,i)./fom_solution(:,i),fespace,'surf',[0 10^(-3)]);
    title('Logarithm of the error');
    
    pause(0.25)
end


