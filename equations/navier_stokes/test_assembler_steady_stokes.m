% common variables
tol = 1e-1;

%% Test 1: verify order convergence against exact solution P2P1

u1ex = @(x) sin(x(2,:)*pi);
u2ex = @(x) exp(x(1,:));
pex = @(x) -0.5*x(1,:)^2;

u1exdx = @(x) 0;
u1exdy = @(x) pi*cos(x(2,:)*pi);
u1exdxdx = @(x) 0;
u1exdydy = @(x) -pi^2*sin(x(2,:)*pi);

u2exdx = @(x) exp(x(1,:));
u2exdy = @(x) 0;
u2exdxdx = @(x) exp(x(1,:));
u2exdydy = @(x) 0;

graduex = @(x) [u1exdx(x) u1exdy(x);u2exdx(x) u2exdy(x)];

pexdx = @(x) -x(1,:);
pexdy = @(x) 0;

mu = 1;
f = @(x) [-mu*(u1exdxdx(x)+u1exdydy(x)) + u1ex(x).*u1exdx(x) + u2ex(x).*u1exdy(x) + pexdx(x);
          -mu*(u2exdxdx(x)+u2exdydy(x)) + u1ex(x).*u2exdx(x) + u2ex(x).*u2exdy(x) + pexdy(x)];

dirichlet_functions = @(x) [u1ex(x) u2ex(x);
                            u1ex(x) u2ex(x);
                            u1ex(x) u2ex(x);
                            u1ex(x) u2ex(x)]';
neumann_functions = @(x) [mu*graduex(x)*[0;-1]-pex(x)*[0;-1], ...
                          mu*graduex(x)*[1;0]-pex(x)*[1;0], ...
                          mu*graduex(x)*[0;1]-pex(x)*[0;1], ...
                          mu*graduex(x)*[-1;0]-pex(x)*[-1;0]];


h = [];
l2errs_u = [];
h1errs_u = [];

l2errs_p = [];

for i = 4:5
    
    n1 = 5*2^(i-1);
    n2 = n1;
    
    h = [h;1/n1];
    
    mesh = create_mesh(0,0,1,1,n1,n2);
    bc_flags = [1 0 1 1];
    fespace_u = create_fespace(mesh,'P2',bc_flags);
    fespace_p = create_fespace(mesh,'P1',bc_flags);
    
    [A,b] = assembler_steady_navier_stokes(fespace_u,fespace_p,f,mu,dirichlet_functions, ...
        neumann_functions);

    % solve stokes problem with u = 0 to get initial guess for newton's method
    u0 = zeros(size(fespace_u.nodes,1)*2,1);
    x0 = A(u0)\b;

    method.name = 'newton';
    method.f = @(u) A(u)*u-b;
    method.x0 = x0;
    method.jac = @(u) build_jac_navier_stokes(A,u,fespace_u);
    method.tol = 1e-8;
    method.maxit = 100;
    
    [sol,err,it] = solve_fluid_system(A,b,fespace_u,fespace_p,method);
    
    l2error_u = compute_L2_error_velocity(fespace_u,sol,@(x) [u1ex(x);u2ex(x)]);
    h1error_u = compute_H1_error_velocity(fespace_u,sol,@(x) [u1ex(x);u2ex(x)], ...
        @(x) [u1exdx(x) u1exdy(x);
        u2exdx(x) u2exdy(x)]);
    
    l2errs_u = [l2errs_u;l2error_u];
    h1errs_u = [h1errs_u;h1error_u];
    
    l2error_p = compute_L2_error(fespace_p,sol.p,pex);
    l2errs_p = [l2errs_p;l2error_p];
    
end

error = h1errs_u + l2errs_p;

c_order = log(error(2:end)./error(1:end-1))./log(h(2:end)./h(1:end-1));

assert(abs(c_order-2) < tol);

%% Test 2: verify order convergence against exact solution P3P2

u1ex = @(x) sin(x(2,:)*pi);
u2ex = @(x) exp(x(1,:));
pex = @(x) -0.5*x(1,:)^2;

u1exdx = @(x) 0;
u1exdy = @(x) pi*cos(x(2,:)*pi);
u1exdxdx = @(x) 0;
u1exdydy = @(x) -pi^2*sin(x(2,:)*pi);

u2exdx = @(x) exp(x(1,:));
u2exdy = @(x) 0;
u2exdxdx = @(x) exp(x(1,:));
u2exdydy = @(x) 0;

graduex = @(x) [u1exdx(x) u1exdy(x);u2exdx(x) u2exdy(x)];

pexdx = @(x) -x(1,:);
pexdy = @(x) 0;

mu = 1;
f = @(x) [-mu*(u1exdxdx(x)+u1exdydy(x)) + u1ex(x).*u1exdx(x) + u2ex(x).*u1exdy(x) + pexdx(x);
          -mu*(u2exdxdx(x)+u2exdydy(x)) + u1ex(x).*u2exdx(x) + u2ex(x).*u2exdy(x) + pexdy(x)];

dirichlet_functions = @(x) [u1ex(x) u2ex(x);
                            u1ex(x) u2ex(x);
                            u1ex(x) u2ex(x);
                            u1ex(x) u2ex(x)]';
neumann_functions = @(x) [mu*graduex(x)*[0;-1]-pex(x)*[0;-1], ...
                          mu*graduex(x)*[1;0]-pex(x)*[1;0], ...
                          mu*graduex(x)*[0;1]-pex(x)*[0;1], ...
                          mu*graduex(x)*[-1;0]-pex(x)*[-1;0]];


h = [];
l2errs_u = [];
h1errs_u = [];

l2errs_p = [];

for i = 4:5
    
    n1 = 5*2^(i-1);
    n2 = n1;
    
    h = [h;1/n1];
    
    mesh = create_mesh(0,0,1,1,n1,n2);
    bc_flags = [1 0 1 1];
    fespace_u = create_fespace(mesh,'P3',bc_flags);
    fespace_p = create_fespace(mesh,'P2',bc_flags);
    
    [A,b] = assembler_steady_navier_stokes(fespace_u,fespace_p,f,mu,dirichlet_functions, ...
        neumann_functions);

    % solve stokes problem with u = 0 to get initial guess for newton's method
    u0 = zeros(size(fespace_u.nodes,1)*2,1);
    x0 = A(u0)\b;

    method.name = 'newton';
    method.f = @(u) A(u)*u-b;
    method.x0 = x0;
    method.jac = @(u) build_jac_navier_stokes(A,u,fespace_u);
    method.tol = 1e-8;
    method.maxit = 100;
    
    [sol,err,it] = solve_fluid_system(A,b,fespace_u,fespace_p,method);
    
    l2error_u = compute_L2_error_velocity(fespace_u,sol,@(x) [u1ex(x);u2ex(x)]);
    h1error_u = compute_H1_error_velocity(fespace_u,sol,@(x) [u1ex(x);u2ex(x)], ...
        @(x) [u1exdx(x) u1exdy(x);
        u2exdx(x) u2exdy(x)]);
    
    l2errs_u = [l2errs_u;l2error_u];
    h1errs_u = [h1errs_u;h1error_u];
    
    l2error_p = compute_L2_error(fespace_p,sol.p,pex);
    l2errs_p = [l2errs_p;l2error_p];
    
end

error = h1errs_u + l2errs_p;

c_order = log(error(2:end)./error(1:end-1))./log(h(2:end)./h(1:end-1));

assert(abs(c_order-3) < tol);