clear all
close all
clc
% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 40;
n1 = 40;

% Create the mesh
mesh = create_mesh(L,H,n1,n2);

% womersley parameters
mu =  1;
rho = 1;
nu = @(x) mu/rho;
f = 14;
w = 2*pi*f;
A = 1;
womersley_number = H/2*sqrt(w*rho/mu);
disp(['Womersley number is = ', num2str(womersley_number)]);

dt_init = 1/(2*f);
dt_init = 1/(20*f);
%dt_init = dt_init/2^(2/7);


% exact solution
rxy   = @(x) abs(x(2)-H/2);
v1_r  = @(r,t) real(A/(1i*rho*w)*(1-besselj(0,2*r/H*womersley_number*1i^(1.5))/(besselj(0,womersley_number*1i^(1.5))))*exp(1i*w*t));
v1    = @(x,t) v1_r(rxy(x),t);
v2    = @(x,t) 0;
p_ex  = @(x,t) cos(w*t)*A/L*(L-x(1));
nullf = @(x) 0;

tt = 0:0.01:10;
x = 0:0.01:H;

f = @(t,x) [0;0];
dirichlet_functions = @(t,x) [0 0;v1(x,t) 0;0 0;v1(x,t) 0]';
neumann_functions = @(t,x) [0 0;0 0;0 0;p_ex(x,t) 0]';

% Create finite element space
bc = [1 0 1 0]; 

poly_degree = 'P2';
fespace_u = create_fespace(mesh,poly_degree,bc);

poly_degree = 'P1';
fespace_p = create_fespace(mesh,poly_degree,bc);

% options
opts.integrate_f = 0;
opts.integrate_neumann = 1;

t0 = 0;
T = 6*dt_init;

n = 3;
total_err = zeros(n,1);

dts = [];

for i = 1:n

    % dt = dt_init/2^((i-1)/(2*n));
    dt = dt_init/2^i
    dts = [dts;dt];
    sol = solver_navier_stokes(fespace_u,fespace_p,t0,T,dt,f,@(x) [v1(x,0);0],@(x) p_ex(x,0),nu,dirichlet_functions,neumann_functions,opts);
    % plot_exact_solution_vs_approx(@(x,t) [v1(x,t);v2(x,t)], @(x,t) p_ex(x,t), sol,'U');

    % check error

    n_timesteps = size(sol.u,2);

    count = 0;
    err = [];
    t = t0;
    while(count < n_timesteps)
        count = count + 1;
        err1 = compute_error(fespace_u,sol.u1(:,count),@(x) v1(x,t),nullf,'L2');
        err2 = compute_error(fespace_u,sol.u2(:,count),@(x) v2(x,t),nullf,'L2');
        err = [err;err1^2+err2^2];
%         err1 = compute_error(fespace_p,sol.p(:,count),@(x) p_ex(x,t),nullf,'L2');
%         err = [err;err1^2];
        t = t + dt;
    end
    
    total_err(i) = sqrt(sum(err)*dt)

end

loglog(dts,total_err,'.-','Linewidth',1,'Markersize',20)
hold on
loglog(dts,dts*1e-2,'--k')
loglog(dts,dts*1e-2*0,'--k')
%axis([min(dts) max(dts) min(total_err) max(total_err)])



