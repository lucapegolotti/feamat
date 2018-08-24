function sol = solve_structure_system(A,b,fespace,varargin)
% Solve a system  of type A sol = b (with A linear matrix) and assemble
% fluid solution with data members
%
% input=
%           A: matrix. It can be constant or anonymous function (A(sol))
%           b: right handside
%           fespace: finite element space of the displacement
%
% output= 
%           sol: data structure with fluid members
%

vecsol = A\b;   

n_nodes = size(fespace.nodes,1);

indices_u1 = 1:n_nodes;
indices_u2 = n_nodes+1:2*n_nodes;

sol.n_nodes = n_nodes;
sol.u1 = vecsol(indices_u1);
sol.u2 = vecsol(indices_u2);
sol.fespace = fespace;
