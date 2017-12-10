function [sol]= solve_fluid_system(A,b,fespace_u,fespace_p)
% Solve a system  of type A sol = b (with A linear matrix) and assemble
% fluid solution with data members
%
% input=
%           A: matrix. It can be constant or anonymous function (A(sol))
%           b: right handside
%           fespace_u: finite element space of the velocity
%           fespace_p: finite element space of the pressure
%
% output= 
%           sol: data structure with fluid members
%

% this is the case when A is a constant matrix
if (~isa(A,'function_handle'))
    vecsol = A\b;
else
    error('The non-linear solver is not yet implemented!');
end

n_nodes_u = size(fespace_u.nodes,1);
n_nodes_p = size(fespace_p.nodes,1);

indices_u1 = 1:n_nodes_u;
indices_u2 = n_nodes_u+1:2*n_nodes_u;
indices_p = 2*n_nodes_u+1:2*n_nodes_u+n_nodes_p;

sol.n_nodes_u = n_nodes_u;
sol.n_nodes_p = n_nodes_p;
sol.u1 = vecsol(indices_u1);
sol.u2 = vecsol(indices_u2);
sol.p = vecsol(indices_p);
sol.fespace_u = fespace_u;
sol.fespace_p = fespace_p;