function [A,b,uLift,INodes] = apply_dirichlet_bc_lifting(A,b,fespace,dirichlet_functions)
% Apply Dirichlet boundary conditions to matrix and rhs using the liftinig
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
%           uLift: lifting Function
%           INodes: internal nodes 

n_nodes = size(A,1);

nodes = fespace.nodes;
bc_flags = fespace.bc;

indices_to_diagonalize = [];
for i = 1:n_nodes
    if (nodes(i,3)~=0 && bc_flags(nodes(i,3)) || ...
        nodes(i,4) ~= 0 && bc_flags(nodes(i,4)))
        indices_to_diagonalize = [indices_to_diagonalize;i];
    end
end

% lifting function
uLift = apply_dirichlet_bc_rhs(zeros(size(b)),fespace,dirichlet_functions);

b = b - A*uLift; 

% eliminate rows and columns
b(indices_to_diagonalize,:) = [];
A(:,indices_to_diagonalize) = [];
A(indices_to_diagonalize,:) = [];

% internal nodes list
INodes = setdiff( [1:n_nodes], indices_to_diagonalize );
