clear all
clc

mesh = read_mesh('bifurcation_separate_boundaries_coarse.msh');

fespace = create_fespace(mesh,'P2',[1;0;1;0;1;1]);

dirichlet_functions = @(x) [0 0; 0 0; 1 1; 0 0; 0 0; 0 0]';
neumann_functions = @(x) [0 0; 0 0; 0 0; 0 0; 0 0; 0 0]';

poisson = 0.2;
young = 10000;

[A,b] = assembler_linear_elasticity(fespace,[0;0],poisson,young, ... 
                                    dirichlet_functions,neumann_functions);
                                
sol = solve_structure_system(A,b,fespace);
[mesh,fespace] = move_mesh(fespace,sol);

draw_mesh(fespace.mesh)

