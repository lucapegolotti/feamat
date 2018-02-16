function [x,err,it] = solve_with_fixed_point(f,x0,tol,maxit)
% Solve non-linear system f(x) = 0 by using fixed point iterations
%
% input=
%           f: anonymous function
%           x0: initial guess for solution
%           tol: tolerance over the error
%           maxit: max number of iterations
%
% output= 
%           sol: solution

res = f(x0);
err = norm(res);

it = 0;
x = x0;
while (err > tol && it < maxit)
    it = it + 1;
    disp(['Fixed point iteration ',num2str(it),' ...']);
    
    x = res + x;
    res = f(x);
    err = norm(res);
    disp(['	done, error = ',num2str(err)]);
end