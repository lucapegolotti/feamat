% clear all
% close all
% clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 40;
n1 = n2*L;

% Create and display the mesh
mesh = create_mesh(L,H,n1,n2);
draw_mesh(mesh);

f = @(x) 1;
mu = @(x) 1;
b = @(x) x*0 + [5;5];
c = @(x) 0;
dirichlet_functions = @(x) [0;0;0;0];
neumann_functions = @(x) [0;0;0;0];

% Create finite element space
bc = [1 1 1 1]; 

poly_degree = 'P2';
fespace = create_fespace(mesh,poly_degree,bc);

% Assemble matrix and rhs
[A,b] =   assembler_diffusion_transport_advection(fespace,f,mu,b,c,dirichlet_functions,neumann_functions);

% Solve the linear system
sol = A\b;

figure
plot_solution_on_fespace(fespace,sol)
pbaspect([1 1 1])

l2norm = compute_norm(fespace,sol,'L2');
display(['Norm = ', num2str(l2norm)]);

%% Here we check the convergence of the error
clear all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

solex = @(x) sin(pi*x(1))*sin(pi*x(2))*x(1)*x(2);
gradex = @(x) [x(2)*sin(pi*x(2))*(sin(pi*x(1)) + pi*x(1)*cos(pi*x(1))); ...
               x(1)*sin(pi*x(1))*(sin(pi*x(2)) + pi*x(2)*cos(pi*x(2)))];

f = @(x) -(2*pi*(x(1)*sin(pi*x(1))*cos(pi*x(2)) + x(2)*sin(pi*x(2))*(cos(pi*x(1)) - pi*x(1)*sin(pi*x(1))))) + ...
           x(2)*sin(pi*x(2))*(sin(pi*x(1)) + pi*x(1)*cos(pi*x(1)));
       

mu = @(x) 1;
b = @(x) x*0 + [1;0];
c = @(x) 0;
dirichlet_functions = @(x) [0;0;0;0];
neumann_functions = @(x) [0;0;0;0];

errl2 = [];
errh1 = [];
h = [];

order = 2;

for i = 1:5

    n2 = 5*2^(i-1);
    n1 = n2*L;

    % Create and display the mesh
    mesh = create_mesh(L,H,n1,n2);
    
    h = [h 1/n2];

    % Create finite element space
    bc = [1 1 1 1]; 

    poly_degree = ['P',num2str(order)];
    fespace = create_fespace(mesh,poly_degree,bc);

    % Assemble matrix and rhs
    [A,rhs] =   assembler_diffusion_transport_advection(fespace,f,mu,b,c,dirichlet_functions,neumann_functions);

    % Solve the linear system
    tic
    disp(['Solution of linear system']);
    sol = A\rhs;
    elapsed = toc;
    disp(['Elapsed time = ', num2str(elapsed),' s']);
    disp('------------------------------');
    
    l2error = compute_error(fespace,sol,solex,gradex,'L2');
    l2norm = compute_norm(fespace,sol,'L2');

    disp(['Relative L2 error = ', num2str(l2error/l2norm)]);
    disp(' ');
    errl2 = [errl2 l2error];
    
    h1error = compute_error(fespace,sol,solex,gradex,'H1');
    h1norm = compute_norm(fespace,sol,'H1');

    disp(['Relative H1 error = ', num2str(h1error/h1norm)]);
    disp(' ');
    errh1 = [errh1 h1error];
    
end

loglog(h,errl2,'.-r','Linewidth',3,'Markersize',20)
hold on
loglog(h,h.^(order+1),'--r','Linewidth',1);

loglog(h,errh1,'.-b','Linewidth',3,'Markersize',20)
hold on
loglog(h,h.^(order),'--b','Linewidth',1);

minl2 = min(errl2);
minh1 = min(errh1);
maxl2 = max(errl2);
maxh1 = max(errh1);

legend('L2 error',['h^',num2str(order+1)],'H1 error',['h^',num2str(order)],'Location','Southeast');
pbaspect([1 1 1]);
axis([min(h) max(h) min([minl2 minh1]) max([maxl2 maxh1])])
set(gca,'Fontsize',25);
