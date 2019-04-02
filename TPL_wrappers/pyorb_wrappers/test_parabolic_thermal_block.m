clear all
close all
clc 

%%

fem_specifics.number_of_elements = 100;
fem_specifics.polynomial_degree = 'P1';
fem_specifics.model = 'thermal_block';
fem_specifics.use_nonhomogeneous_dirichlet = 'N';
fem_specifics.mesh_name = 'cube100x100';

% New fields in time dependent case!!!!
fem_specifics.number_of_time_instances = 200;
fem_specifics.final_time = 1.0;
fem_specifics.theta = 1.0;

params = [1.00, 1.00, 1.00]; 

sol = solve_parameter( params, fem_specifics );


%%
[~, fespace] = set_fem_simulation( fem_specifics );
x = fespace.nodes(:,1);
y = fespace.nodes(:,2);
T = fem_specifics.final_time;
exact_sol = exp(-T) * sin(2*pi*x) .* sin(2*pi*y);

%%
err = norm(exact_sol - sol.u(:,end),2);

%%
figure
plot3(x,y,exact_sol);
grid on
title('Exact solution');

figure
plot3(x,y,sol.u(:,end));
grid on
title('Numerical solution');

%%
figure
plot3(x,y,exact_sol - sol.u(:,end));
grid on
title('Pointwise error');