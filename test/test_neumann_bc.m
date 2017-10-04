clear allc
clc

% dummy test on small vector

mesh = create_mesh(0,0,1,1,2,2);

bc = [0 0 0 0];

fespace = create_fespace(mesh,'P1',bc);
 
neumann_bc =@(x) [1;0;1;0];

b = zeros(size(fespace.nodes,1),1);

b = apply_neumann_bc(fespace,b,neumann_bc);

%% actual test on poisson problem

L = 1;
H = 1;

n2 = 100;
n1 = n2*L;

f = @(x) 0;
mu = @(x) 1;
dirichlet_functions = @(x) [0;0;0;x(2)*(1-x(2))];
neumann_functions = @(x) [0;(x(2)-1)*x(2);0;0];

% Create finite element space
bc = [0 0 0 1]; 

poly_degree = 'P2';
fespace = create_fespace(mesh,poly_degree,bc);

% Assemble matrix and rhs
[A,b] =   assembler_poisson(fespace,f,mu,dirichlet_functions,neumann_functions);

% Solve the linear system
sol = A\b;

h = 0.001;
xp1 = L-h;
xp2 = L;
[x,y1] = get_values_over_line(fespace,sol,100,xp1,'Ypar');
[~,y2] = get_values_over_line(fespace,sol,100,xp2,'Ypar');

%plot_solution_on_fespace(fespace,sol)

plot(x,(y2-y1)/h)
hold on
plot(x,(x-1).*x)



