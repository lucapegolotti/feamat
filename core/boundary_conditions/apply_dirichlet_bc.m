function [A,b] = apply_dirichlet_bc(A,b,fespace,dirichlet_functions)
% Apply Dirichlet boundary conditions to matrix and rhs 
%
% input=
%           A: matrix
%           b: right handside
%           fespace: finite element space
%           dirichlet_functions: anonymous function describing the boundary
%                                data
% output=
%           A: diagonalized matrix
%           b: right handside with boundary data

A = apply_dirichlet_bc_matrix(A,fespace,1);
b = apply_dirichlet_bc_rhs(b,fespace,dirichlet_functions);