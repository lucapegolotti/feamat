function [x,err,it] = solve_with_newtons_method(f,x0,jac,tol,maxit)
% Solve non-linear system f(x) = 0 by using Newton's method
%
% input=
%           f: anonymous function
%           x0: initial guess for solution
%           jac: jacobian of the function
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
    disp(['Newton iteration ',num2str(it),' ...']);
    J = jac(x);
    x = x - J\res;
    res = f(x);
    err = norm(res);
    disp(['	done, error = ',num2str(err)]);
end



