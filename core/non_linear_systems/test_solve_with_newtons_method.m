% common variables
tol = 1e-12;

%% Test 1: test with scalar equations

f = @(x) sin(x);
df = @(x) cos(x);

x0 = 2.5;

x = solve_with_newtons_method(f,x0,df,tol*0.1,100);

assert(abs(x-pi) < tol)

%% Test 2: test with system of equations

f_aux = @(x) [sin(x(1)*x(2));exp(x(1))];
f = @(x) f_aux(x) - f_aux([3;2]);
jac = @(x) [cos(x(1)*x(2))*x(2) cos(x(1)*x(2))*x(1); exp(x(1)) 0];
  
x0 = [4;1.5];

[x,err,it] = solve_with_newtons_method(f,x0,jac,tol*0.1,100);

assert(norm(x-[3;2]) < tol)

