clear all
close all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 30;
n1 = 30;

% Create and display the mesh
mesh = create_mesh(L,H,n1,n2);

f = @(x) [0;0];
nu = @(x) 1;
dirichlet_functions = @(x) [0 0;0 0;0 0;0 1]';
neumann_functions = @(x) [0 0;0 0;0 0;0 0]';

% Create finite element space
bc = [1 0 1 1]; 

poly_degree = 'P2';
fespace_u = create_fespace(mesh,poly_degree,bc);

poly_degree = 'P1';
fespace_p = create_fespace(mesh,poly_degree,bc);

% Solve problem
[sol] = solver_steady_navier_stokes(fespace_u,fespace_p,f,nu,dirichlet_functions,neumann_functions);

subplot(1,2,1)
plot_solution_vp(fespace_u,fespace_p,sol.u,'U');
title('V')
axis([0 L 0 H])

pbaspect([L H 1])
subplot(1,2,2)
plot_solution_vp(fespace_u,fespace_p,sol.u,'P')
title('P')
pbaspect([L H 1])

%% Here we check the convergence of the error
clear all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

solexu1 = @(x) sin(x(2)*pi);
gradexu1 = @(x) [0; pi*cos(x(2)*pi)];

solexu2 = @(x) 0;
gradexu2 = @(x) [0;0];

solexp = @(x) -0.5*x(1)^2+0.5;

f = @(x) [pi^2*sin(x(2)*pi)-x(1);0];
nu = @(x) 1;
dirichlet_functions = @(x) [0 0;0 0;0 0;solexu1(x) 0]';
neumann_functions = @(x) [0 0;0 0;0 0;0 0]';

errul2 = [];
erruh1 = [];
errpl2 = [];

h = [];

for i = 1:5

    n2 = 5*2^(i-1);
    n1 = n2*L;

    % Create and display the mesh
    mesh = create_mesh(L,H,n1,n2);
    
    h = [h 1/n2];

    % Create finite element space
    bc = [1 0 1 1]; 

    fespace_u = create_fespace(mesh,'P2',bc);
    fespace_p = create_fespace(mesh,'P1',bc);
    
    n_nodes_u = size(fespace_u.nodes,1);
    n_nodes_p = size(fespace_p.nodes,1);

    % Solve problem
    [sol] = solver_steady_navier_stokes(fespace_u,fespace_p,f,nu,dirichlet_functions,neumann_functions);

    u1 = sol.u1;
    u2 = sol.u2;
    p = sol.p;
    elapsed = toc;
    disp(['Elapsed time = ', num2str(elapsed),' s']);
    disp('------------------------------');
    
    err1 = compute_error(fespace_u,u1,solexu1,gradexu1,'L2');
    err2 = compute_error(fespace_u,u2,solexu2,gradexu2,'L2');
    
    l2error = sqrt(err1^2 + err2^2);
    norm1 = compute_norm(fespace_u,u1,'L2');
    norm2 = compute_norm(fespace_u,u2,'L2');
    l2norm = sqrt(norm1^2 + norm2^2);
    
    disp(['Relative L2 error = ', num2str(l2error/l2norm)]);
    disp(' ');
    errul2 = [errul2 l2error];
    
    err1 = compute_error(fespace_u,u1,solexu1,gradexu1,'H1');
    err2 = compute_error(fespace_u,u2,solexu2,gradexu2,'H1');
    
    h1error = sqrt(err1^2 + err2^2);
    norm1 = compute_norm(fespace_u,u1,'H1');
    norm2 = compute_norm(fespace_u,u2,'H1');
    h1norm = sqrt(norm1^2 + norm2^2);

    disp(['Relative H1 error = ', num2str(h1error/h1norm)]);
    disp(' ');
    erruh1 = [erruh1 h1error];
    
    err1 = compute_error(fespace_p,p,solexp,@(x) x*0,'L2');
    pnorm = compute_norm(fespace_p,p,'L2');

    disp(['Relative L2 error on pressure = ', num2str(err1/pnorm)]);
    disp(' ');
    errpl2 = [errpl2 err1];
    
end

subplot(1,2,1)

loglog(h,errul2,'.-r','Linewidth',3,'Markersize',20)
hold on
loglog(h,h.^3,'--r','Linewidth',1);

loglog(h,erruh1,'.-b','Linewidth',3,'Markersize',20)
hold on
loglog(h,h.^2,'--b','Linewidth',1);

minl2 = min(errul2);
minh1 = min(erruh1);
maxl2 = max(errul2);
maxh1 = max(erruh1);

legend('L2 error',['h^',num2str(3)],'H1 error',['h^',num2str(2)],'Location','Southeast');
pbaspect([1 1 1]);
axis([min(h) max(h) min([minl2 minh1]) max([maxl2 maxh1])])
set(gca,'Fontsize',25);

subplot(1,2,2)

loglog(h,errpl2,'.-r','Linewidth',3,'Markersize',20)
hold on
loglog(h,h.^2,'--r','Linewidth',1);

minl2 = min(errpl2);
maxl2 = max(errpl2);

legend('L2 error',['h^',num2str(2)],'Location','Southeast');
pbaspect([1 1 1]);
axis([min(h) max(h) min(minl2,h(end)^2) max(maxl2,h(1)^2)])
set(gca,'Fontsize',25);

