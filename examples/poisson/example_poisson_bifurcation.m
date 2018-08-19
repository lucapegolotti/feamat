clear all
close all
clc

mesh = read_mesh('bifurcation_long.msh');

bc_flags = [0 1 1];
fespace = create_fespace(mesh,'P2',bc_flags);

f = @(x) 1;
mu = 3.5;

dirichlet_functions = @(x) [1;1;1;1]*sin(x(1)*x(2));
neumann_functions = @(x) [0;0;0;0;0];

[A,b] = assembler_poisson(fespace,f,mu,dirichlet_functions,neumann_functions);

sol = A\b;

plot_fe_function(sol,fespace);
axis([mesh.xp mesh.xp+mesh.L mesh.yp mesh.yp+mesh.H])
