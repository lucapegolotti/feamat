%% check convergence in time of rosenbrock methods
clear all
close all
clc

co = [0 0 1;
    0 0.5 0;
    1 0 0;
    0 0.75 0.75;
    0.75 0 0.75;
    0.75 0.75 0;
    0.25 0.25 0.25];
set(groot,'defaultAxesColorOrder',co)

nmethods = 1;
methods = {};
%methods{1} = 'BDF1';
methods{1} = 'ROSI2P2';
methods{2} = 'ROS3Pw';
methods{3} = 'ROS2';
methods{4} = 'ROS3P';
methods{5} = 'ROWDAIND2';
methods{6} = 'ROS34PW2';
methods{7} = 'ROS3PRw';

npoints = 4;
total_errs_ul2 = zeros(npoints,nmethods);
total_errs_pl2 = zeros(npoints,nmethods);
total_errs_uh1 = zeros(npoints,nmethods);


% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 10;
n1 = 10;

% Create and display the mesh
mesh = create_mesh(L,H,n1,n2);

exact_solution = 3;

if (exact_solution == 1)
% Exact solution
uex = @(t,x) sin(t*pi)*(x(2)-H).^3;
vex = @(t,x) t*cos(x(1)*pi);

% derivatives of exact solution
udx = @(t,x) 0;
udy = @(t,x) 3*sin(t*pi)*(x(2)-H).^2;
vdx = @(t,x) -pi*t*sin(x(1)*pi);
vdy = @(t,x) 0;
udxdt = @(t,x) 0;
udydt = @(t,x) 3*pi*cos(t*pi)*(x(2)-H).^2;
vdxdt = @(t,x) -pi^2*t*cos(x(1)*pi);
vdydt = @(t,x) 0;

udxdx = @(t,x) 0;
udxdxdt = @(t,x) 0;
udydy = @(t,x) 6*sin(t*pi)*(x(2)-H);
udydydt = @(t,x) pi*6*cos(t*pi)*(x(2)-H);

vdxdx = @(t,x) -pi^2*t*cos(x(1)*pi);
vdxdxdt = @(t,x) -pi^2*cos(x(2)*pi);
vdydy = @(t,x) 0;
vdydydt = @(t,x) 0;

udt = @(t,x) pi*cos(t*pi)*(x(2)-H).^3;
udtdt = @(t,x) -pi^2*sin(t*pi)*(x(2)-H).^3;
vdt = @(t,x) cos(x(1)*pi);
vdtdt = @(t,x) 0;

w = 2;
pex = @(t,x) (L-x(1))/L*cos(t*w);
pdt = @(t,x) -w*(L-x(1))/L*sin(t*w);
pdx = @(t,x) -cos(t*w)/L;
pdy = @(t,x) 0;
pdxdt = @(t,x) w*sin(t*w)/L;
pdydt = @(t,x) 0;

f = @(t,x) [udt(t,x) - udxdx(t,x) - udydy(t,x) + uex(t,x)*udx(t,x) + ...
    vex(t,x)*udy(t,x) + pdx(t,x);
    vdt(t,x) - vdxdx(t,x) - vdydy(t,x) + uex(t,x)*vdx(t,x) + ...
    vex(t,x)*vdy(t,x) + pdy(t,x)];
f_dt = @(t,x) [udtdt(t,x) - udxdxdt(t,x) - udydydt(t,x) + udt(t,x)*udx(t,x) + ...
    uex(t,x)*udxdt(t,x) + vdt(t,x)*udy(t,x) * vex(t,x)*udydt(t,x) + pdxdt(t,x);
    vdtdt(t,x) - vdxdxdt(t,x) - vdydydt(t,x) + udt(t,x)*vdx(t,x) + ...
    uex(t,x)*vdxdt(t,x) + vdt(t,x)*vdy(t,x) + vex(t,x)*vdydt(t,x) + pdydt(t,x)];

nu = @(x) 1;

dirichlet_functions = @(t,x) [uex(t,x) vex(t,x);
    uex(t,x) vex(t,x);
    uex(t,x) vex(t,x);
    uex(t,x) vex(t,x)]';
dirichlet_dt = @(t,x) [udt(t,x) vdt(t,x);
    udt(t,x) vdt(t,x);
    udt(t,x) vdt(t,x);
    udt(t,x) vdt(t,x)]';

neumann_functions = @(t,x) [-udy(t,x) -vdy(t,x)+pex(t,x);
                            udx(t,x)-pex(t,x) vdx(t,x);
                            udy(t,x) vdy(t,x)-pex(t,x);
                            -udx(t,x)+pex(t,x) -vdx(t,x)]';
neumann_dt = @(t,x) [-udydt(t,x) -vdydt(t,x)+pdt(t,x);
                     udxdt(t,x)-pdt(t,x) vdxdt(t,x);
                     udydt(t,x) vdydt(t,x)-pdt(t,x);
                     -udxdt(t,x)+pdt(t,x) -vdxdt(t,x)]';
elseif (exact_solution == 2)
% Exact solution
uex = @(t,x) sin(t)*(x(2)^2+x(1));
vex = @(t,x) sin(t)*(x(1)^2-x(2));

% derivatives of exact solution
udx = @(t,x) sin(t);
udy = @(t,x) sin(t)*(2*x(2));
vdx = @(t,x) sin(t)*(2*x(1));
vdy = @(t,x) -sin(t);
udxdt = @(t,x) cos(t);
udydt = @(t,x) cos(t)*(2*x(2));
vdxdt = @(t,x) cos(t)*(2*x(1));
vdydt = @(t,x) -cos(t);

udxdx = @(t,x) 0;
udxdxdt = @(t,x) 0;
udydy = @(t,x) 2*sin(t);
udydydt = @(t,x) 2*cos(t);

vdxdx = @(t,x) 2*sin(t);
vdxdxdt = @(t,x) 2*cos(t);
vdydy = @(t,x) 0;
vdydydt = @(t,x) 0;

udt = @(t,x) cos(t)*(x(2)^2+x(1));
udtdt = @(t,x) -sin(t)*(x(2)^2+x(1));
vdt = @(t,x) cos(t)*(x(1)^2-x(2));
vdtdt = @(t,x) -sin(t)*(x(1)^2-x(2));

w = 2;
pex = @(t,x) exp(-t)*(x(1)+x(2)-1);
pdt = @(t,x) -exp(-t)*(x(1)+x(2)-1);
pdx = @(t,x) exp(-t);
pdy = @(t,x) exp(-t);
pdxdt = @(t,x) -exp(-t);
pdydt = @(t,x) -exp(-t);

f = @(t,x) [udt(t,x) - udxdx(t,x) - udydy(t,x) + uex(t,x)*udx(t,x) + ...
    vex(t,x)*udy(t,x) + pdx(t,x);
    vdt(t,x) - vdxdx(t,x) - vdydy(t,x) + uex(t,x)*vdx(t,x) + ...
    vex(t,x)*vdy(t,x) + pdy(t,x)];
f_dt = @(t,x) [udtdt(t,x) - udxdxdt(t,x) - udydydt(t,x) + udt(t,x)*udx(t,x) + ...
    uex(t,x)*udxdt(t,x) + vdt(t,x)*udy(t,x) * vex(t,x)*udydt(t,x) + pdxdt(t,x);
    vdtdt(t,x) - vdxdxdt(t,x) - vdydydt(t,x) + udt(t,x)*vdx(t,x) + ...
    uex(t,x)*vdxdt(t,x) + vdt(t,x)*vdy(t,x) + vex(t,x)*vdydt(t,x) + pdydt(t,x)];

nu = @(x) 1;

dirichlet_functions = @(t,x) [uex(t,x) vex(t,x);
    uex(t,x) vex(t,x);
    uex(t,x) vex(t,x);
    uex(t,x) vex(t,x)]';
dirichlet_dt = @(t,x) [udt(t,x) vdt(t,x);
    udt(t,x) vdt(t,x);
    udt(t,x) vdt(t,x);
    udt(t,x) vdt(t,x)]';

neumann_functions = @(t,x) [-udy(t,x) -vdy(t,x)+pex(t,x);
    udx(t,x)-pex(t,x) vdx(t,x);
    udy(t,x) vdy(t,x)-pex(t,x);
    -udx(t,x)+pex(t,x) -vdx(t,x)]';
neumann_dt = @(t,x) [-udydt(t,x) -vdydt(t,x)+pdt(t,x);
    udxdt(t,x)-pdt(t,x) vdxdt(t,x);
    udydt(t,x) vdydt(t,x)-pdt(t,x);
    -udxdt(t,x)+pdt(t,x) -vdxdt(t,x)]';
elseif (exact_solution == 3)
% Exact solution
uex = @(t,x) sin(t)*(x(2)*(H-x(2)));
vex = @(t,x) 0;

% derivatives of exact solution
udx = @(t,x) 0;
udy = @(t,x) sin(t)*(H-2*x(2));
vdx = @(t,x) 0;
vdy = @(t,x) 0;
udxdt = @(t,x) 0;
udydt = @(t,x) cos(t)*(H-2*x(2));
vdxdt = @(t,x) 0;
vdydt = @(t,x) 0;

udxdx = @(t,x) 0;
udxdxdt = @(t,x) 0;
udydy = @(t,x) -2*sin(t);
udydydt = @(t,x) -2*cos(t);

vdxdx = @(t,x) 0;
vdxdxdt = @(t,x) 0;
vdydy = @(t,x) 0;
vdydydt = @(t,x) 0;

udt = @(t,x) cos(t)*(x(2)*(H-x(2)));
udtdt = @(t,x) -sin(t)*(x(2)*(H-x(2)));
vdt = @(t,x) 0;
vdtdt = @(t,x) 0;

w = 2*pi;
pex = @(t,x) cos(w*t)*(L-x(1))/L;
pdt = @(t,x) -w*sin(w*t)*(L-x(1))/L;
pdx = @(t,x) -cos(w*t)/L;
pdy = @(t,x) 0;
pdxdt = @(t,x) -w*sin(w*t)/L;
pdydt = @(t,x) 0;

f = @(t,x) [udt(t,x) - udxdx(t,x) - udydy(t,x) + uex(t,x)*udx(t,x) + ...
    vex(t,x)*udy(t,x) + pdx(t,x);
    vdt(t,x) - vdxdx(t,x) - vdydy(t,x) + uex(t,x)*vdx(t,x) + ...
    vex(t,x)*vdy(t,x) + pdy(t,x)];
f_dt = @(t,x) [udtdt(t,x) - udxdxdt(t,x) - udydydt(t,x) + udt(t,x)*udx(t,x) + ...
    uex(t,x)*udxdt(t,x) + vdt(t,x)*udy(t,x) * vex(t,x)*udydt(t,x) + pdxdt(t,x);
    vdtdt(t,x) - vdxdxdt(t,x) - vdydydt(t,x) + udt(t,x)*vdx(t,x) + ...
    uex(t,x)*vdxdt(t,x) + vdt(t,x)*vdy(t,x) + vex(t,x)*vdydt(t,x) + pdydt(t,x)];

nu = @(x) 1;

dirichlet_functions = @(t,x) [uex(t,x) vex(t,x);
    uex(t,x) vex(t,x);
    uex(t,x) vex(t,x);
    uex(t,x) vex(t,x)]';
dirichlet_dt = @(t,x) [udt(t,x) vdt(t,x);
    udt(t,x) vdt(t,x);
    udt(t,x) vdt(t,x);
    udt(t,x) vdt(t,x)]';

neumann_functions = @(t,x) [-udy(t,x) -vdy(t,x)+pex(t,x);
    udx(t,x)-pex(t,x) vdx(t,x);
    udy(t,x) vdy(t,x)-pex(t,x);
    -udx(t,x)+pex(t,x) -vdx(t,x)]';
neumann_dt = @(t,x) [-udydt(t,x) -vdydt(t,x)+pdt(t,x);
    udxdt(t,x)-pdt(t,x) vdxdt(t,x);
    udydt(t,x) vdydt(t,x)-pdt(t,x);
    -udxdt(t,x)+pdt(t,x) -vdxdt(t,x)]';
end

% Create finite element space
bc = [1 1 1 1];

poly_degree = 'P2';
fespace_u = create_fespace(mesh,poly_degree,bc);

poly_degree = 'P1';
fespace_p = create_fespace(mesh,poly_degree,bc);

% options
opts.integrate_f = 1;
opts.integrate_neumann = 1;


dt_init = 0.05;
t0 = 0;

final_time_error = 1;

for method = 1:nmethods
    total_err_ul2 = zeros(npoints,1);
    total_err_uh1 = zeros(npoints,1);
    total_err_pl2 = zeros(npoints,1);
    dts = [];
    for i = 1:npoints
        
        dt = dt_init/2^((i-1));
        dts = [dts;dt];
        T = dt_init*2;
        
        if (strcmp(methods{method},'BDF1'))
            sol = solver_navier_stokes(fespace_u,fespace_p,t0,T,dt,f,@(x) [uex(t0,x);vex(t0,x)],@(x) pex(t0,x),nu,dirichlet_functions,neumann_functions,opts);
        else
            sol = solver_navier_stokes_rosenbrock(fespace_u,fespace_p,t0,T,dt,f, ...
                f_dt,@(x) [uex(t0,x);vex(t0,x)],@(x) pex(t0,x),nu,dirichlet_functions, ...
                dirichlet_dt,neumann_functions,neumann_dt,methods{method},opts);
        end
        
        
        if (final_time_error)
            % check error
            err1 = compute_error(fespace_u,sol.u1(:,end),@(x) uex(T,x),@(x) [udx(T,x);udy(T,x)],'L2');
            err2 = compute_error(fespace_u,sol.u2(:,end),@(x) vex(T,x),@(x) [vdx(T,x);vdy(T,x)],'L2');
            
            total_err_ul2(i) = sqrt(err1^2+err2^2);
          
            
            err1 = compute_error(fespace_u,sol.u1(:,end),@(x) uex(T,x),@(x) [udx(T,x);udy(T,x)],'H1');
            err2 = compute_error(fespace_u,sol.u2(:,end),@(x) vex(T,x),@(x) [vdx(T,x);vdy(T,x)],'H1');
            
            total_err_uh1(i) = sqrt(err1^2+err2^2);
            
            diff = sol.p(1,end) - pex(T,[0;0]);
            
            total_err_pl2(i) = compute_error(fespace_p,sol.p(:,end)-diff,@(x) pex(T,x),@(x) 0,'L2');
        else
            ntimesteps = size(sol.u,2);
            
            count = 0;
            
            t = t0;
            errul2 = 0;
            erruh1 = 0;
            errpl2 = 0;
            
            while (count < ntimesteps)
                count = count + 1;
                err1 = compute_error(fespace_u,sol.u1(:,count),@(x) uex(t,x),@(x) [udx(t,x);udy(t,x)],'L2');
                err2 = compute_error(fespace_u,sol.u2(:,count),@(x) vex(t,x),@(x) [vdx(t,x);vdy(t,x)],'L2');
                
                errul2 = errul2 + err1^2+err2^2;
                
                err1 = compute_error(fespace_u,sol.u1(:,count),@(x) uex(t,x),@(x) [udx(t,x);udy(t,x)],'H1');
                err2 = compute_error(fespace_u,sol.u2(:,count),@(x) vex(t,x),@(x) [vdx(t,x);vdy(t,x)],'H1');
                
                erruh1 = erruh1 + err1^2+err2^2;
                
                errpl2 = errpl2 + compute_error(fespace_p,sol.p(:,count),@(x) pex(t,x),@(x) 0,'L2')^2;
                t = t + dt;
            end
            total_err_ul2(i) = sqrt(dt*errul2);
            total_err_uh1(i) = sqrt(dt*erruh1);
            total_err_pl2(i) = sqrt(dt*errpl2);
        end
    end
    total_errs_ul2(:,method) = total_err_ul2;
    total_errs_uh1(:,method) = total_err_uh1;
    total_errs_pl2(:,method) = total_err_pl2;
    
end

%%
clc
close all
%load data
% plot l2 errors
methods{nmethods+1}='\Delta t';
methods{nmethods+2}='\Delta t^2';
methods{nmethods+3}='\Delta t^3';

figure(1)
loglog(dts,total_errs_ul2,'.-')
title('L2 error')

hold on

loglog(dts,dts,'--k');
loglog(dts,dts.^2,'-.k');
loglog(dts,dts.^3*1,':k');

log(total_errs_ul2(2:end,:)./total_errs_ul2(1:end-1,:))./log(dts(2:end)./dts(1:end-1))

legend(methods)

% plot h1 errors

figure(2)
loglog(dts,total_errs_uh1,'.-')
title('H1 error')

hold on

loglog(dts,dts,'--k');
loglog(dts,dts.^2,'-.k');
loglog(dts,dts.^3*1,':k');

log(total_errs_uh1(2:end,:)./total_errs_uh1(1:end-1,:))./log(dts(2:end)./dts(1:end-1))

legend(methods)

% plot l2 error on p

figure(3)
loglog(dts,total_errs_pl2,'.-')
title('l2 error on p')

hold on

loglog(dts,dts,'--k');
loglog(dts,dts.^2,'-.k');
loglog(dts,dts.^3*1,':k');

log(total_errs_pl2(2:end,:)./total_errs_pl2(1:end-1,:))./log(dts(2:end)./dts(1:end-1))

legend(methods)

%%
plot_exact_solution_vs_approx(@(x,t) [uex(t,x);vex(t,x)],@(t,x) pex(x,t), sol, 'U_dif')