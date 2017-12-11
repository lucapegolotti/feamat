function err = compute_H1_error_velocity(fespace,sol,fexact,gradexact)
% Compute H1 error of fluid solution with respect to exact solution of the
% velocity
% input=
%           fespace: finite element space
%           sol: fluid data structure
%           fexact: exact solution
% output=
%           err: L2 error
%

err1 = compute_H1_error(fespace,sol.u1,@(x) fexact(x)'*[1;0], ...
                        @(x) gradexact(x)'*[1;0]);
err2 = compute_H1_error(fespace,sol.u2,@(x) fexact(x)'*[0;1], ...
                        @(x) gradexact(x)'*[0;1]);

err = sqrt(err1^2+err2^2);