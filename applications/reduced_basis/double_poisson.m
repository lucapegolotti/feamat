%% solving the poisson problem on two separate domains with rb
%  The solution is glued togheter by imposing the continuity of the
%  integrals

clear all
clc

% create mesh on left subdomain
xp1 = 0;
yp1 = 0;
L1 = 0.5;
H1 = 1;

mesh1 = create_mesh(xp1,yp1,L1,H1,20,30);

% create mesh on right subdomain
xp2 = 0.5;
yp2 = 0;
L2 = 0.5;
H2 = 1;

mesh2 = create_mesh(xp2,yp2,L2,H2,20,20);

% draw global mesh
meshes = {};
meshes{end+1} = mesh1;
meshes{end+1} = mesh2;
%draw_multimesh(meshes);

fun = @(x) sin(x(1))*cos(x(2));

% coefficients for parameters
rangemus = [0.1 10];
nmus = 5;
rangeomegas = [-10 10];
nomegas = 5;

% take snapshots in left subdomain
fespace1 = create_fespace(mesh1,'P1',[1 0 1 1]);
[Sn1,A1,b1,v1,v1t] = take_snapshots_poisson(fespace1,fun,rangemus,nmus,rangeomegas,nomegas,2);
% plot_solution_on_fespace(fespace1,Sn1(:,5))
% hold on

% take snapshots in right subdomain
fespace2 = create_fespace(mesh2,'P1',[1 1 1 0]);
[Sn2,A2,b2,v2,v2t] = take_snapshots_poisson(fespace2,fun,rangemus,nmus,rangeomegas,nomegas,4);
% plot_solution_on_fespace(fespace2,Sn2(:,2))
% axis([0 1 0 1])

% compute svd 
[V1,S1,~] = svd(Sn1);
[V2,S2,~] = svd(Sn2);

% select modes such that 99.9% energy is preserved
d1 = diag(S1);
totalsum = sum(d1);
partialsum = 0;

count = 0;
while (partialsum < 0.999 * totalsum)
    count = count + 1;
    partialsum = partialsum + d1(count);
end

V1 = V1(:,1:count);

disp(['Selecting ',num2str(count),' modes in left subdomain. Energy = ', num2str(partialsum/totalsum*100), '%']);

d2 = diag(S2);
totalsum = sum(d2);
partialsum = 0;

count = 0;
while (partialsum < 0.999 * totalsum)
    count = count + 1;
    partialsum = partialsum + d2(count);
end

V2 = V2(:,1:count);

disp(['Selecting ',num2str(count),' modes in right subdomain. Energy = ', num2str(partialsum/totalsum*100), '%']);

% solving for value of parameter mu = pi
mu = pi;

n1 = size(A1,1);
n2 = size(A2,2);

Ah = [A1 zeros(n1,n2) v1;
      zeros(n2,n1) A2 -v2;
      v1t -v2t 0];
  
nmodes1 = size(V1,2);
nmodes2 = size(V2,2);
nnodes1 = size(V1,1);
nnodes2 = size(V2,1);
W = [V1 zeros(nnodes1,nmodes2) zeros(nnodes1,1);
     zeros(nnodes2,nmodes1) V2 zeros(nnodes2,1);
     zeros(1,nmodes1) zeros(1,nmodes2) 1];
 
Arb = W'*Ah*W;

rhs = [V1'*b1;V2'*b2;0];

totalsol = (mu*Arb)\rhs;

sol1 = totalsol(1:nmodes1);
sol2 = totalsol(nmodes1+1:nmodes1+nmodes2);

plot_solution_on_fespace(fespace1,V1*sol1)
hold on

plot_solution_on_fespace(fespace2,V2*sol2)

% compute exact solution

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 30;
n1 = 30;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);

mu = @(x) pi;
dirichlet_functions = @(x) [0;0;0;0];
neumann_functions = @(x) [0;0;0;0];

% Create finite element space
bc = [1 1 1 1]; 

poly_degree = 'P1';
fespace = create_fespace(mesh,poly_degree,bc);

% Assemble matrix and rhs
[A,b] =   assembler_poisson(fespace,fun,mu,dirichlet_functions,neumann_functions);

% Solve the linear system
sol = A\b;

figure
plot_solution_on_fespace(fespace,sol)
pbaspect([1 1 1])



