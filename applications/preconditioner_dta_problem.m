clear all
close all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 10;
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

bc2 = [0 0 0 0]; 

poly_degree = 'P2';
fespace2 = create_fespace(mesh,poly_degree,bc2);

% Assemble matrix and rhs
[M,~] =   assembler_diffusion_transport_advection(fespace2,f,mu,@(x) b(x)*0,@(x) 1,dirichlet_functions,neumann_functions);
[A,b] =   assembler_diffusion_transport_advection(fespace,f,mu,b,c,dirichlet_functions,neumann_functions);

%M = eye(size(A));

% M = assemble_mass(fespace);
% M = apply_dirichlet_bc_matrix(M,fespace,1);
% M = apply_dirichlet_bc_matrix(M',fespace,1);
% 
% Solve the linear system
tic
[L,U] = ilu(A);
%sol = pcg(A,b,1e-8,100,L,U);
sol = gmres(A,b,100,1e-8,100,L,U);
toc

% modified system
tic
H = ichol(M);
Pm = H'*H;
Pa = L*U;
Pls = Pa'*(M\(Pa));
aux1 = Pm\A;
mat = A'*aux1;
%cond(A)
%cond(Pls\mat)
%cond(mat)
aux2 = Pm\b;
rhs = A'*aux2;
% sol = gmres(mat,rhs,100,1e-8,1000,Pls);
% sol = gmres(mat,rhs,100,1e-8,1000
sol = gmres(mat,rhs,1000000,1e-7,1000,Pls);
toc

figure
plot_solution_on_fespace(fespace,sol)
pbaspect([1 1 1])

l2norm = compute_norm(fespace,sol,'L2');
display(['Norm = ', num2str(l2norm)]);

( sqrt(cond(Pls\mat)) - 1 ) / ( sqrt(cond(Pls\mat)) + 1 )
( sqrt(cond(mat)) - 1 ) / ( sqrt(cond(mat)) + 1 )
( sqrt(cond(Pa\mat)) - 1 ) / ( sqrt(cond(Pa\mat)) + 1 )

