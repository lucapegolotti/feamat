function [A,b] = apply_dirichlet_bc(A,b,fespace,dirichlet_functions, varargin)
% Apply Dirichlet boundary conditions to matrix and rhs 
%
% input=
%           A: matrix
%           b: right handside
%           fespace: finite element space
%           dirichlet_functions: anonymous function describing the boundary
%                                data
%           varargin: contains time, if unsteady problem is solved
% output=
%           A: diagonalized matrix
%           b: right handside with boundary data

A = apply_dirichlet_bc_matrix(A,fespace,1);
if nargin < 5
    b = apply_dirichlet_bc_rhs(b,fespace,dirichlet_functions);
else
    b = apply_dirichlet_bc_rhs(b,fespace,dirichlet_functions, varargin{1});
end