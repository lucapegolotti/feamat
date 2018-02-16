function [sol,err,it]= solve_fluid_system(A,b,fespace_u,fespace_p,varargin)
% Solve a system  of type A sol = b (with A linear matrix) and assemble
% fluid solution with data members
%
% input=
%           A: matrix. It can be constant or anonymous function (A(sol))
%           b: right handside
%           fespace_u: finite element space of the velocity
%           fespace_p: finite element space of the pressure
%           (required if A is function)
%           'method' structure containing method's parameters
%
% output= 
%           sol: data structure with fluid members
%           err: error (for non-linear systems)
%           it: number of iterations (for iterative solver)
%

% this is the case when A is a constant matrix
if (~isa(A,'function_handle'))
    vecsol = A\b;
    err = 0;
    it  = 0;
else
    if (nargin < 5)
        error('The method to use should be specified as optional parameter!');
    end
    method = varargin{1};
    if (strcmp(method.name,'newton'))
        [vecsol,err,it] = solve_with_newtons_method(method.f,method.x0,method.jac,method.tol,method.maxit);
    elseif (strcmp(method.name,'fixed_point'))
        [vecsol,err,it] = solve_with_fixed_point(method.f,method.x0,method.tol,method.maxit);
    elseif (strcmp(method.name,'newton_and_fixed_point'))
        [vecsol,err,it] = solve_with_newtons_method(method.f,method.x0,method.jac,method.tol,method.maxit);
        if (err > tol)
            oldvecsol = vecsol;
            [vecsol,err,it] = solve_with_newtons_method(method.f,vecsol,method.jac,method.tol,method.maxit);
            if (it > method.maxit)
                vecsol = oldvecsol;
                warning('Neither Newton method nor fixed point converged with the desired tolerance');
            end
        end
    else
        error('The non-linear solver is not yet implemented!');
    end
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