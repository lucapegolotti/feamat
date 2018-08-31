function [A] = apply_dirichlet_bc_matrix(A,fespace,diagvalue)
% Apply Dirichlet boundary conditions to matrix  putting to zero
% the rows of the matrix corresponding to Dirichlet nodes and replacing the
% diagonal values by diagvalue
% input=
%           A: matrix
%           fespace: finite element space
%           diagvalue: diagonal value to insert in the matrix
%
% output=
%           A: diagonalized matrix

indices_to_diagonalize = find_dirichlet_indices(fespace);

n_indices = length(indices_to_diagonalize);

A(indices_to_diagonalize,:) = 0;
A(indices_to_diagonalize,indices_to_diagonalize) = spdiags(ones(n_indices,1)*diagvalue, ...
                                                   0,n_indices,n_indices);


