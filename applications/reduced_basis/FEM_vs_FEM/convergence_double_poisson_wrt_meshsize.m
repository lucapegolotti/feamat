%clear all
%clc

errs = [];

for N = [20 40 80 160]
% solving the poisson problem on two subdomains with FEM
%  The final solution is glued by imposing the continuity

fun = @(x) sin(x(1))*cos(10*x(2))*exp(x(1))+log(x(2))*x(1);
mu = pi;

% L2 error of full solution
disp('Computing exact solution');

mesh = create_mesh(0,0,1,1,N,N);

fespace = create_fespace(mesh,'P2',[1 1 1 1]);
[A,b] = assembler_poisson(fespace,fun,@(x) mu,@(x) [0;0;0;0],@(x) [0;0;0;0]);

totalsol = A\b;


% create mesh on left subdomain
xp1 = 0;
yp1 = 0;
L1 = 0.5;
H1 = 1;

n1x = N/2;
n2x = N/2;
n1y = N;
n2y = N;

mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);

% create mesh on right subdomain
xp2 = 0.5;
yp2 = 0;
L2 = 0.5;
H2 = 1;

mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);

% draw global mesh
meshes = {};
meshes{end+1} = mesh1;
meshes{end+1} = mesh2;
%draw_multimesh(meshes);

% solve problems on subdomains
mu = pi;

fespace1 = create_fespace(mesh1,'P2',[1 0 1 1]);
fespace2 = create_fespace(mesh2,'P2',[1 1 1 0]);

[A1,b1] = assembler_poisson(fespace1,fun,@(x) mu,@(x) [0;0;0;0],@(x) [0;0;0;0]);
[A2,b2] = assembler_poisson(fespace2,fun,@(x) mu,@(x) [0;0;0;0],@(x) [0;0;0;0]);

n1 = size(A1,1);
n2 = size(A2,1);

indices1 = 1:n1;
indices2 = n1+1:n1+n2;

V1 = [];
V2 = [];

nmodes = 16;

l2err = [];

for i = 1:nmodes
    disp(['Solving with mode with frequency omega * ',num2str(i)]);
    v1 = zeros(n1,1);
    v2 = zeros(n2,1);
    
    v1 = apply_neumann_bc(fespace1,v1,@(x) [0;sin(x(2)*pi*i);0;0]);
    v2 = apply_neumann_bc(fespace2,v2,@(x) [0;0;0;sin(x(2)*pi*i)]);
    v11 = apply_neumann_bc(fespace1,v1,@(x) [0;cos(x(2)*pi*i);0;0]);
    v22 = apply_neumann_bc(fespace2,v2,@(x) [0;0;0;cos(x(2)*pi*i)]);
    
    V1 = [V1 v1 v11];
    V2 = [V2 v2 v22];
    
    mat = [A1 sparse(n1,n2) -V1; sparse(n2,n1) A2 V2; -V1' V2' sparse(i*2,i*2)];
    mat = sparse(mat);
    
    mat(indices1,:) = apply_dirichlet_bc_matrix(mat(indices1,:),fespace1,1);
    mat(indices2,n1+1:end) = apply_dirichlet_bc_matrix(mat(indices2,n1+1:end),fespace2,1);
    
    rhs = [b1;b2;zeros(2*i,1)];
    
    sol = mat\rhs;
    
    sol1 = sol(indices1);
    sol2 = sol(indices2);
    
    figure(1);
    plot_solution_on_fespace(fespace1,sol1)
    hold on

    plot_solution_on_fespace(fespace2,sol2)
    
    hold off
    
    interpsol = interp_2solutions_on_fespace(fespace1,sol1,fespace2,sol2,fespace);
    err = compute_error(fespace,abs(totalsol-interpsol),@(x)0,@(x)[0;0],'L2');
    disp(['Total L2 error = ', num2str(err)]);
    l2err = [l2err;err];
    % check convergence of derivatives
    yy = 0:0.01:1;
    x1 = 0.499;
    x2 = 0.5;
    
    dudx = zeros(length(yy),1);
    
    figure(2)
    for j = 1:length(dudx)
        dudx(j) = (interpolate_in_point(fespace1,sol1,x2,yy(j))-interpolate_in_point(fespace1,sol1,x1,yy(j)))/(x2-x1);
    end
    plot(yy,dudx)
    hold on
    dudx = zeros(length(yy),1);
    x1 = 0.5;
    x2 = 0.501;
    
    for j = 1:length(dudx)
        dudx(j) = (interpolate_in_point(fespace2,sol2,x2,yy(j))-interpolate_in_point(fespace2,sol2,x1,yy(j)))/(x2-x1);
    end
    plot(yy,dudx)
    
    
    dudx = zeros(length(yy),1);
    x1 = 0.4995;
    x2 = 0.5005;
    for j = 1:length(dudx)
        dudx(j) = (interpolate_in_point(fespace,totalsol,x2,yy(j))-interpolate_in_point(fespace,totalsol,x1,yy(j)))/(x2-x1);
    end
    plot(yy,dudx)
    title('du/dx over the interface')
    hold off
    
    % check convergence of function at interface
    x1 = 0.5;
    
    u = zeros(length(yy),1);
    
    figure(3)
    for j = 1:length(u)
        u(j) = interpolate_in_point(fespace1,sol1,x1,yy(j));
    end
    plot(yy,u)
    hold on
    
    u = zeros(length(yy),1);
    
    for j = 1:length(u)
        u(j) = interpolate_in_point(fespace2,sol2,x1,yy(j));
    end
    plot(yy,u)
    

    u = zeros(length(yy),1);

    for j = 1:length(u)
        u(j) = interpolate_in_point(fespace,totalsol,x1,yy(j));
    end
    plot(yy,u)
    title('u at the interface')
    hold off

    pause(0.1)
end

errs = [errs l2err];

figure(4)
semilogy(1:nmodes,l2err,'.-','Markersize',10)
title('Convergence l2 error')
xlabel('number of modes')
ylabel('l2 error')

end

save('errors_convergence_wrt_h.mat',errs)