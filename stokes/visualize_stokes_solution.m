function visualize_stokes_solution(solution, delay,varargin)
print = 0;
if (nargin >= 3)
    print = 1;
    namefile = varargin{1};
end

u = solution.u;
t0 = solution.t0;
T = solution.T;
dt = solution.dt;
mesh = solution.mesh;
fespace_u = solution.fespace_u;
fespace_p = solution.fespace_p;
L = mesh.L;
H = mesh.H;
n_timesteps = size(u,2);
t = t0;
count = 0;
figure(1)
pause
while (count < n_timesteps)  
    count = count + 1;
    subplot(1,2,1)
    plot_solution_vp(fespace_u,fespace_p,u(:,count),'U',[solution.minnorm solution.maxnorm])
    title(['V at t = ',num2str(t)])
    axis([0 L 0 H])

    pbaspect([L H 1])
    subplot(1,2,2)
    plot_solution_vp(fespace_u,fespace_p,u(:,count),'P',[solution.minp solution.maxp])
    title(['P at t = ',num2str(t)])
    pbaspect([L H 1])
    pause(delay)
    if (print)
        saveas(gcf,['data/',namefile','_',num2str(count),'.png'])
    end
    t = t + dt;
end