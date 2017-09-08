function [solution] = solver_navier_stokes(fespace_u,fespace_p,t0,T,dt,fun,u0,p0,nu,dirichlet_functions,neumann_functions,varargin)
integrate_f = 1;
if (nargin >= 11)
    opts = varargin{1};
    integrate_f = opts.integrate_f;
    integrate_neumann = opts.integrate_neumann;
end


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

zero_mat_u = zeros(n_nodes_u);
zero_mat_p = zeros(n_nodes_p);
zero_mat_up = zeros(n_nodes_u,n_nodes_p);

M = [m zero_mat_u zero_mat_up; 
     zero_mat_u m zero_mat_up; 
     zero_mat_up' zero_mat_up' zero_mat_p];

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
    t = t + dt;

    disp(['Time = ',num2str(t)]);

    fun1 = @(x) fun(t,x)'*[1;0];
    fun2 = @(x) fun(t,x)'*[0;1];
    dir1 = @(x) dirichlet_functions(t,x)'*[1;0];
    dir2 = @(x) dirichlet_functions(t,x)'*[0;1];
    neu1 = @(x) neumann_functions(t,x)'*[1;0];
    neu2 = @(x) neumann_functions(t,x)'*[0;1];

    % C = assemble_convective_term(fespace_u,u);
    % C = apply_dirichlet_bc_matrix(C,fespace_u,0);
    
%     C = @(u) apply_dirichlet_bc_matrix(assemble_convective_term(fespace_u,u), ... 
%                                        fespace_u,0);
    
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
    
    C = assemble_convective_term(fespace_u,[u(nodes_u1);u(nodes_u2)])*1000;
    % apply_dirichlet_bc_matrix(C,fespace_u,0);
    
    H1 = [A+C zero_mat_u -B1_u];
    H2 = [zero_mat_u A+C -B2_u];
    H3 = [-B1 -B2 zero_mat_p];
    H = [H1;H2;H3];

    mat = 1/dt * M + H;
    mat(nodes_u1,nodes_u1) = apply_dirichlet_bc_matrix(mat(nodes_u1,nodes_u1),fespace_u,1);
    mat(nodes_u2,nodes_u2) = apply_dirichlet_bc_matrix(mat(nodes_u2,nodes_u2),fespace_u,1);

    b = [b1;b2;zeros(n_nodes_p,1)];
    
    rhs = b + 1/dt * M * u;
    
%     b1 = apply_dirichlet_bc_rhs(b1,fespace_u,dir1);
%     b2 = apply_dirichlet_bc_rhs(b2,fespace_u,dir2);
    
    rhs(nodes_u1) = apply_dirichlet_bc_rhs(rhs(nodes_u1),fespace_u,dir1);
    rhs(nodes_u2) = apply_dirichlet_bc_rhs(rhs(nodes_u2),fespace_u,dir2);
    % u = fixed_point(u,n_nodes_u,fespace_u,M,H,b,dt,1e-9,10);
    
    u = mat\rhs;
    
    u1 = u(1:n_vertices);
    u2 = u(n_nodes_u+1:n_nodes_u+1+n_vertices);
    
    minnorm = min(min(sqrt(u1(1:n_vertices).^2+u2(1:n_vertices).^2)),minnorm);
    maxnorm = max(max(sqrt(u1(1:n_vertices).^2+u2(1:n_vertices).^2)),maxnorm);
    minp = min(min(u(2*n_nodes_u+1:end)),minp);
    maxp = max(max(u(2*n_nodes_u+1:end)),maxp);
    
    sol(:,count) = u;
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
