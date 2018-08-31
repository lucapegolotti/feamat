
%   Author: Stefano Pagani <stefano.pagani at polimi.it>

clear all
clc

% geometry  definition
bottom_left_corner_x = 0;
bottom_left_corner_y = 0;

L = 1.0;
H = 1.0;

% number of elements
n_elements_x = 40;
n_elements_y = 40;

mesh = create_mesh(bottom_left_corner_x, ...
                   bottom_left_corner_y, ...
                   L,H,n_elements_x,n_elements_y);

% boundary conditions                
bc_flags = [1 1 1 1];
dirichlet_functions = @(x) [0,0,0,0];
neumann_functions = @(x) [0;0;0;0];

% finite elements space definition
fespace   = create_fespace(mesh,'P2',bc_flags);
fespaceP1 = create_fespace(mesh,'P1',bc_flags);

% forcing term
f = @(x) 1.0 + 0.*x(1,:).*x(2,:);

% covariance definition
std_a = 0.1;
const = 1;

numRV = 8;

% compute correlation matrix
for k_1=1:length(mesh.vertices)
    for k_2=k_1:length(mesh.vertices)
        C(k_1,k_2) = std_a^2 * exp( - abs( mesh.vertices(k_1,1) - mesh.vertices(k_2,1) )/const - abs( mesh.vertices(k_1,2) - mesh.vertices(k_2,2) )/const );
        C(k_2,k_1) = C(k_1,k_2);
    end
end

% compute eigen functions
[V,D] = eigs(C./(std_a^2), [], numRV);

KL_EXP     = V*sqrt(D(1:numRV,1:numRV));



% conductivity mean-value 
mu = 1.0;

% Assembler
[A_mean,b, uL, iN] = assembler_poisson_lifting(fespace,f,mu,dirichlet_functions,neumann_functions);
        

for i_KL=1:numRV
    
    % K-L mode 
    mu = @(x) evaluate_fe_function(std_a*KL_EXP(:,i_KL),fespaceP1,x);

    % Assembly
    [A_in{i_KL},b_in{i_KL}] = assembler_poisson_lifting(fespace,0,mu,dirichlet_functions,neumann_functions);
    

end



% random sample
psi = randn(numRV,1);


% linear combination
A = A_mean;    
for i_KL=1:numRV        
    A = A + psi(i_KL)*A_in{i_KL};
    b = b + psi(i_KL)*b_in{i_KL};
end


% Solver
sol = uL;
sol(iN)  = A\b;

% plot of the solution
plot_fe_function(sol,fespace)
%export_vtk_scalar(sol,fespace,'example_thermal_block.vtk');