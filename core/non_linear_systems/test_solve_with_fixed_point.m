% common variables
tol = 1e-12;

%% Test 1: test with scalar equations

f = @(x) sin(x);

x0 = 3;

x = solve_with_fixed_point(f,x0,tol*0.1,100);

assert(abs(x-pi) < tol)
