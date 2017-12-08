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

for i = 1:n_nodes
    if (nodes(i,3)~=0)
        if (bc_flags(nodes(i,3)))
            A(i,:) = 0;
            if (i <= dim2)
                A(i,i) = diagvalue;
            end
        elseif (nodes(i,4) ~= 0 && bc_flags(nodes(i,4)))
            A(i,:) = 0;
            if (i <= dim2)
                A(i,i) = diagvalue;
            end
        end
    end
end
