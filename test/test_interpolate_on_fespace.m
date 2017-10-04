clear all
close all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 30;
n1 = n2*L;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);
draw_mesh(mesh);

f = @(x) 0;
mu = @(x) 1;
dirichlet_functions = @(x) [0;0;0;x(2)*(1-x(2))];
neumann_functions = @(x) [0;0;0;0];

% Create finite element space
bc = [1 0 1 1]; 

poly_degree = 'P1';
fespace = create_fespace(mesh,poly_degree,bc);

% Assemble matrix and rhs
[A,b] =   assembler_poisson(fespace,f,mu,dirichlet_functions,neumann_functions);

% Solve the linear system
sol = A\b;

interpsol = interp_on_fespace(fespace,sol,fespace);

norm(sol-sol)

%% Convergence with respect to numerical solution
clear all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

f = @(x) x(2)^2*log(x(1));
mu = @(x) 1;
dirichlet_functions = @(x) [0;0;0;0];
neumann_functions = @(x) [-pi*sin(pi*x(1)).*cos(pi*x(2));
                          pi*cos(pi*x(1)).*sin(pi*x(2));
                          pi*sin(pi*x(1))*cos(pi*x(2));
                          -pi*cos(pi*x(1)).*sin(pi*x(2))];
                      
% create reference solution

bc = [1 0 1 0]; 

order = 2;

mesh = create_mesh(0,0,1,1,140,140);
fespace_fine = create_fespace(mesh,['P',num2str(order)],bc);
[A,b] = assembler_poisson(fespace_fine,f,mu,dirichlet_functions,neumann_functions);

globalsol = A\b;
l2norm = compute_norm(fespace_fine,globalsol,'L2');


errl2 = [];
h = [];

for i = 1:5

    n2 = 2*2^(i-1);
    n1 = n2*L;

    % Create and display the mesh
    mesh = create_mesh(0,0,L,H,n1,n2);
    
    h = [h 1/n2];

    poly_degree = ['P',num2str(order)];
    fespace = create_fespace(mesh,poly_degree,bc);

    % Assemble matrix and rhs
    [A,b] = assembler_poisson(fespace,f,mu,dirichlet_functions,neumann_functions);

    % Solve the linear system
    tic
    disp(['Solution of linear system']);
    sol = A\b;
    elapsed = toc;
    disp(['Elapsed time = ', num2str(elapsed),' s']);
    disp('------------------------------');
    
    interpsol = interp_on_fespace(fespace,sol,fespace_fine);
    
    l2error = compute_error(fespace_fine,(globalsol-interpsol),@(x) 0,@(x) [0;0],'L2');

    disp(['Relative L2 error = ', num2str(l2error/l2norm)]);
    disp(' ');
    errl2 = [errl2 l2error];
    
end

loglog(h,errl2,'.-r','Linewidth',3,'Markersize',20)
hold on
loglog(h,h.^(order+1),'--r','Linewidth',1);

minl2 = min(errl2);
maxl2 = max(errl2);

legend('L2 error',['h^',num2str(order+1)],'Location','Southeast');
pbaspect([1 1 1]);
axis([min(h) max(h) minl2 maxl2])
set(gca,'Fontsize',25);