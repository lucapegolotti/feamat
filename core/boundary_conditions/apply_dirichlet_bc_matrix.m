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


n_nodes = size(A,1);
dim2 = size(A,2);

nodes = fespace.nodes;
bc_flags = fespace.bc;

indices_to_diagonalize = [];
% for i = 1:n_nodes
%     if (nodes(i,3)~=0)
%         if (bc_flags(nodes(i,3)))
%             A(i,:) = 0;
%             if (i <= dim2)
%                 A(i,i) = diagvalue;
%             end
%         elseif (nodes(i,4) ~= 0 && bc_flags(nodes(i,4)))
%             A(i,:) = 0;
%             if (i <= dim2)
%                 A(i,i) = diagvalue;
%             end
%         end
%     end
% end
indices_to_diagonalize = [];
for i = 1:n_nodes
    if (nodes(i,3)~=0 && bc_flags(nodes(i,3)) || ...
        nodes(i,4) ~= 0 && bc_flags(nodes(i,4)))
        indices_to_diagonalize = [indices_to_diagonalize;i];
    end
end

n_indices = length(indices_to_diagonalize);

A(indices_to_diagonalize,:) = 0;
A(indices_to_diagonalize,indices_to_diagonalize) = spdiags(ones(n_indices,1)*diagvalue, ...
                                                   0,n_indices,n_indices);


