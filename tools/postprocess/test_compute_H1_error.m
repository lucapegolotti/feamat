% common variables 
tol = 1e-12;

%% Test 1: check if optimized code gives same result

solex = @(x) sin(pi*x(1)).*sin(pi*x(2));
gradex = @(x) [cos(pi*x(1)).*sin(pi*x(2));cos(pi*x(2)).*sin(pi*x(1))]*pi;

mesh = create_mesh(0,0,1,1,20,20);
fespace = create_fespace(mesh,'P2',[1 1 1 1]);

f = @(x) 2*pi^2*sin(pi*x(1,:)).*sin(pi*x(2,:));
mu = 1;
dirichlet_functions = @(x) [0;0;0;0];
neumann_functions = @(x) [-pi*sin(pi*x(1)).*cos(pi*x(2));
    pi*cos(pi*x(1)).*sin(pi*x(2));
    pi*sin(pi*x(1))*cos(pi*x(2));
    -pi*cos(pi*x(1)).*sin(pi*x(2))];

[A,b] = assembler_poisson(fespace,f,mu,dirichlet_functions,neumann_functions);

sol = A\b;

h1error1 = compute_H1_error(fespace,sol,solex,gradex);

fespace.mesh.type = '';
h1error2 = compute_H1_error(fespace,sol,solex,gradex);

assert(abs((h1error1 - h1error2)) < tol)