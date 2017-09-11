function [solution] = solver_navier_stokes_rosenbrock(fespace_u,fespace_p,t0,T,dt,fun,fundt,u0,p0,nu,dirichlet_functions, dirichlet_dt, neumann_functions, neumann_dt, method,varargin)
integrate_f = 1;
integrate_neumann = 1;

if (nargin >= 16)
    opts = varargin{1};
    integrate_f = opts.integrate_f;
    integrate_neumann = opts.integrate_neumann;
end

% load rosenbrock coefficients
coeffs = rosenbrock_coeffs(method);
nstages = coeffs.nstages;

bc_flags_u = fespace_u.bc;
thereisneumann = 1;

n_nodes_u = size(fespace_u.nodes,1);
n_nodes_p = size(fespace_p.nodes,1);

nodes_u1 = 1:n_nodes_u;
nodes_u2 = n_nodes_u+1:n_nodes_u*2;

if (length(find(bc_flags_u)) == 4)
    thereisneumann = 0;
end

% set initial condition
u1 = project_function(fespace_u,@(x) u0(x)'*[1;0]);
u2 = project_function(fespace_u,@(x) u0(x)'*[0;1]);
p  = project_function(fespace_p,@(x) p0(x));

u = [u1;u2;p];

% assemble constant matrices
A = assemble_stiffness(nu,fespace_u);
B1 = assemble_divergence(fespace_u,fespace_p,'dx');
B2 = assemble_divergence(fespace_u,fespace_p,'dy');
m = assemble_mass(fespace_u);

A = apply_dirichlet_bc_matrix(A,fespace_u,1);

zero_mat_u = zeros(n_nodes_u);
zero_mat_p = zeros(n_nodes_p);
zero_mat_up = zeros(n_nodes_u,n_nodes_p);

M = [m zero_mat_u zero_mat_up;
    zero_mat_u m zero_mat_up;
    zero_mat_up' zero_mat_up' zero_mat_p];
M = sparse(M);

% A = apply_dirichlet_bc_matrix(A,fespace_u,1);
B1_u = B1';
B2_u = B2';
B1_u = apply_dirichlet_bc_matrix(B1_u,fespace_u,0);
B2_u = apply_dirichlet_bc_matrix(B2_u,fespace_u,0);

t = t0;
sol = zeros(n_nodes_u*2+n_nodes_p,ceil((T-t0)/dt));
sol(:,1) = u;

n_vertices = size(fespace_u.mesh.vertices,1);

count = 1;
minnorm = min(sqrt(u1(1:n_vertices).^2+u2(1:n_vertices).^2));
maxnorm = max(sqrt(u1(1:n_vertices).^2+u2(1:n_vertices).^2));
minp = inf;
maxp = 0;
while (T-t>dt/2)
    count = count + 1;
    
    disp(['Time = ',num2str(t+dt)]);
    
    C = assemble_convective_term(fespace_u,u);
    J11 = assemble_Jac_convective_term(fespace_u,u,1,1);
    J12 = assemble_Jac_convective_term(fespace_u,u,1,2);
    J21 = assemble_Jac_convective_term(fespace_u,u,2,1);
    J22 = assemble_Jac_convective_term(fespace_u,u,2,2);
    
    C = apply_dirichlet_bc_matrix(C,fespace_u,0);
    J11 = apply_dirichlet_bc_matrix(J11,fespace_u,0);
    J12 = apply_dirichlet_bc_matrix(J12,fespace_u,0);
    J21 = apply_dirichlet_bc_matrix(J21,fespace_u,0);
    J22 = apply_dirichlet_bc_matrix(J22,fespace_u,0);

    H1 = [A+C+J11 J12 -B1_u];
    H2 = [J21 A+C+J22 -B2_u];
    H3 = [-B1 -B2 zero_mat_p];
    J = -[H1;H2;H3];
    J = sparse(J);
        
    mat = M - coeffs.gamma(1,1)*dt*J;
    mat(nodes_u1,nodes_u1) = apply_dirichlet_bc_matrix(mat(nodes_u1,nodes_u1),fespace_u,1);
    mat(nodes_u2,nodes_u2) = apply_dirichlet_bc_matrix(mat(nodes_u2,nodes_u2),fespace_u,1);
        
    stages = zeros(length(u),nstages);
    
    for j = 1:nstages
        
        alphai = sum(coeffs.alpha(j,:));
        gammai = sum(coeffs.gamma(j,:));
        
        uhat = u + dt*stages(:,1:j-1)*coeffs.alpha(j,1:j-1)';
        sumstages = stages(:,1:j-1)*coeffs.gamma(j,1:j-1)';
        
        % handle rhs
        fun1 = @(x) (fun(t+alphai*dt,x) + gammai*dt*fundt(t,x))'*[1;0];
        fun2 = @(x) (fun(t+alphai*dt,x) + gammai*dt*fundt(t,x))'*[0;1];
        neu1 = @(x) (neumann_functions(t+alphai*dt,x) + gammai*dt*neumann_dt(t,x))'*[1;0];
        neu2 = @(x) (neumann_functions(t+alphai*dt,x) + gammai*dt*neumann_dt(t,x))'*[0;1];
        dir1 = @(x) (dirichlet_functions(t+alphai*dt,x) + gammai*dt*dirichlet_dt(t,x))'*[1;0];
        dir2 = @(x) (dirichlet_functions(t+alphai*dt,x) + gammai*dt*dirichlet_dt(t,x))'*[0;1];

        Chat = assemble_convective_term(fespace_u,uhat);
        Chat = apply_dirichlet_bc_matrix(Chat,fespace_u,0);
        
        b1 = zeros(n_nodes_u,1);
        b2 = zeros(n_nodes_u,1);
        
        if (integrate_f)
            b1 = assemble_rhs(fespace_u,fun1);
            b2 = assemble_rhs(fespace_u,fun2);
        end
        
        if (thereisneumann && integrate_neumann)
            b1 = apply_neumann_bc(fespace_u,b1,neu1);
            b2 = apply_neumann_bc(fespace_u,b2,neu2);
        end

        matrhs = [A+Chat zero_mat_u -B1_u; 
                  zero_mat_u A+Chat -B2_u;
                  -B1 -B2 zero_mat_p];

        matrhs = sparse(matrhs);

        b = [b1;b2;zeros(n_nodes_p,1)];

        b = b - matrhs*uhat + dt*J*sumstages;

        dirb1 = zeros(length(b1),1);
        dirb2 = zeros(length(b2),1);
        index_dir = zeros(length(b1),1);
        
        dirb1 = apply_dirichlet_bc_rhs(dirb1,fespace_u,dir1);
        dirb2 = apply_dirichlet_bc_rhs(dirb2,fespace_u,dir2);
        index_dir = apply_dirichlet_bc_rhs(index_dir,fespace_u,@(x) [1;1;1;1]);

        for k = 1:n_nodes_u
            if (index_dir(k) == 1)
               b(k) = (dirb1(k)  - uhat(k) - dt*sumstages(k))/(dt*coeffs.gamma(1,1));
               b(k + n_nodes_u) = (dirb2(k) - uhat(k + n_nodes_u) - dt*sumstages(k + n_nodes_u))/(dt*coeffs.gamma(1,1));
            end
        end
        stages(:,j) = mat\b;
    end
    
    u = u + dt * stages * coeffs.b;
    
    % post process
    u1 = u(1:n_vertices);
    u2 = u(n_nodes_u+1:n_nodes_u+1+n_vertices);
    
    minnorm = min(min(sqrt(u1(1:n_vertices).^2+u2(1:n_vertices).^2)),minnorm);
    maxnorm = max(max(sqrt(u1(1:n_vertices).^2+u2(1:n_vertices).^2)),maxnorm);
    minp = min(min(u(2*n_nodes_u+1:end)),minp);
    maxp = max(max(u(2*n_nodes_u+1:end)),maxp);
    
    sol(:,count) = u;
    t = t + dt;
end

solution.u = sol;
solution.u1 = sol(1:n_nodes_u,:);
solution.u2 = sol(n_nodes_u+1:n_nodes_u*2,:);
solution.p = sol(n_nodes_u*2+1:end,:);
solution.t0 = t0;
solution.T = T;
solution.dt = dt;
solution.mesh = fespace_u.mesh;
solution.fespace_u = fespace_u;
solution.fespace_p = fespace_p;
solution.minnorm = minnorm;
solution.maxnorm = maxnorm;
solution.minp = minp;
solution.maxp = maxp;
