function [A] = apply_dirichlet_bc_matrix_indices(A,diagvalue,indices_to_diagonalize)
% If i is in indices_to_diagonalize, this function sets A(i,i) = diagvalue
% and A(i,j) = 0 is i is different from j.
% input=
%           A: matrix
%           diagvalue: diagonal value to insert in the matrix
%           indices_to_diagonalize: list of indices to diagonalize
% output=
%           A: diagonalized matrix

n_indices = length(indices_to_diagonalize);

A(indices_to_diagonalize,:) = 0;

if (diagvalue ~= 0)
    A(indices_to_diagonalize,indices_to_diagonalize) = ...
        spdiags(ones(n_indices,1)*diagvalue,0,n_indices,n_indices);
end