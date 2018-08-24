% common variables
tol = 5e-2;

%% Test 1: verify order convergence against exact solution P2P1
poisson = 0.1;
young = 100;

lambda = poisson * young/((1 - 2*poisson)*(1+poisson));
mu = young/(2*(1+poisson));

u1ex = @(x) sin(x(1,:).*x(2,:));
u2ex = @(x) x(1,:).^2;

u1exdx = @(x) x(2,:).*cos(x(1,:).*x(2,:));
u1exdy = @(x) x(1,:).*cos(x(1,:).*x(2,:));
u1exdxdx = @(x) -x(2,:).^2.*sin(x(1,:).*x(2,:));
u1exdydy = @(x) -x(1,:).^2.*sin(x(1,:).*x(2,:));
u1exdxdy = @(x) -x(1,:).*x(2,:).*sin(x(1,:).*x(2,:)) + cos(x(1,:).*x(2,:));

u2exdx = @(x) 2*x(1,:);
u2exdy = @(x) 0;
u2exdxdx = @(x) 2*x(1,:).^0;
u2exdydy = @(x) 0;
u2exdxdy = @(x) 0;

stress = @(x) [lambda * (u1exdx(x) + u2exdy(x)) + 2*mu*u1exdxdx(x)  mu*(u1exdy(x) + u2exdxdy(x));
                mu*(u1exdxdy(x) + u2exdxdx(x)) lambda * (u1exdx(x) + u2exdy(x)) + 2*mu*u2exdydy(x)];

f = @(x) -[lambda * (u1exdxdx(x) + u2exdxdy(x)) + 2*mu*u1exdxdx(x) + (u1exdydy(x) + u2exdxdy(x)) * mu;
          mu * (u1exdxdy(x) + u2exdxdx(x)) + lambda*(u1exdxdy(x) + u2exdydy(x)) + 2*mu*u2exdydy(x)];

dirichlet_functions = @(x) [u1ex(x) u2ex(x);
    u1ex(x) u2ex(x);
    u1ex(x) u2ex(x);
    u1ex(x) u2ex(x)]';
neumann_functions = @(x) [stress(x)*[0;-1], ...
    stress(x)*[1;0], ...
    stress(x)*[0;1], ...
    stress(x)*[-1;0]];


h = [];
h1errs = [];

for i = 4:5
    
    n1 = 5*2^(i-1);
    n2 = n1;
    
    h = [h;1/n1];
    
    mesh = create_mesh(0,0,1,1,n1,n2);
    bc_flags = [1 1 1 1];
    fespace = create_fespace(mesh,'P2',bc_flags);
    
    [A,b] = assembler_linear_elasticity(fespace,f,poisson,young, ...
        dirichlet_functions, neumann_functions);
    sol = solve_structure_system(A,b,fespace);
    
    h1error = compute_H1_error_velocity(fespace,sol,@(x) [u1ex(x);u2ex(x)], ...
        @(x) [u1exdx(x) u1exdy(x);
        u2exdx(x) u2exdy(x)]);
       
    h1errs = [h1errs;h1error];
end

error = h1errs;

c_order = log(error(2:end)./error(1:end-1))./log(h(2:end)./h(1:end-1));

assert(abs(c_order-2) < tol);