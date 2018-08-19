clear all
close all
clc

tol = 1e-2;

boundary_indicators = @(x) [(norm(x) < 0.3 + tol);
                            abs(x(1)-0.5) < tol;
                            abs(x(1)+0.5) < tol;
                            abs(x(2)-0.5) < tol;
                            abs(x(2)+0.5) < tol];

mesh = read_mesh('square_hole_fine.msh',boundary_indicators);

bc_flags = [1 1 1 1 1];
fespace = create_fespace(mesh,'P2',bc_flags);

f = @(x) sin(x(1,:).*x(2,:)).*x(2,:).^3;
mu = 3.5;

dirichlet_functions = @(x) [sin(x(1)*5)*0.1;-x(2).^2+0.5^2;0;0;0];
neumann_functions = @(x) [0;0;0;0;0];

[A,b] = assembler_poisson(fespace,f,mu,dirichlet_functions,neumann_functions);

sol = A\b;

plot_fe_function(sol,fespace);
axis square
