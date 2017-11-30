function plot_steady_solution_navier_stokes(solution,label)

u = solution.u;

mesh = solution.mesh;
fespace_u = solution.fespace_u;
fespace_p = solution.fespace_p;
vertices = fespace_u.mesh.vertices;
L = mesh.L;
H = mesh.H;

figure(1)
pause
while (count < n_timesteps)  
    count = count + 1;
    subplot(1,2,1)
    if (label == 'U')
        plot_solution_vp(fespace_u,fespace_p,u(:,count),'U',[solution.minnorm solution.maxnorm])
        title(['V at t = ',num2str(t)])
    elseif(label == 'P')
        plot_solution_vp(fespace_u,fespace_p,u(:,count),'P',[solution.minp solution.maxp])
        title(['P at t = ',num2str(t)])
    else
        error('Label is not valid! Must be either U or P');
    end
    axis([0 L 0 H])
    pbaspect([L H 1])
    
    subplot(1,2,2)

    n_vertices = size(fespace_u.mesh.vertices,1);
    n1 = size(fespace_u.mesh.X,1);
    n2 = size(fespace_u.mesh.X,2);
    
    if (label == 'U')
        u1 = zeros(n_vertices,1);
        u2 = zeros(n_vertices,1);

        for i = 1:n_vertices
            u1(i) = abs(vexact(vertices(i,1:2),t)'*[1;0]-u(i,count));
            u2(i) = abs(vexact(vertices(i,1:2),t)'*[0;1]-u(i+size(fespace_u.nodes,1),count));
        end

        U1 = reshape(u1,n1,n2);
        U2 = reshape(u2,n1,n2);

        N = sqrt(U1.^2+U2.^2);

        [~,c] = contourf(fespace_u.mesh.X,fespace_p.mesh.Y,N);
        title(['V exact at t = ',num2str(t)])
        c.LineStyle = 'none';
        shading interp
        h = colorbar;
        m = solution.minnorm;
        M = solution.maxnorm;
        set(h,'YLim',[m*0.9 M*1.1]);
        caxis([ m*0.9 M*1.1 ]);
        c.LevelList = linspace(m,M,20);
        hold on

        q = quiver(fespace_u.mesh.X,fespace_u.mesh.Y,U1,U2);
        q.Color = 'black';
        q.LineWidth = 1;
        q.AutoScaleFactor = 0.8;
        hold off
        axis([0 L 0 H])
        pbaspect([L H 1])
    elseif (label == 'P')
        p = zeros(n_vertices,1);

        for i = 1:n_vertices
            p(i) = pexact(vertices(i,1:2),t);
        end

        P = reshape(p,n1,n2);

        [~,c] = contourf(fespace_u.mesh.X,fespace_p.mesh.Y,P);
        title(['P exact at t = ',num2str(t)])
        c.LineStyle = 'none';
        shading interp
        h = colorbar;
        m = solution.minp;
        M = solution.maxp;
        set(h,'YLim',[m*0.9 M*1.1]);
        caxis([ m*0.9 M*1.1 ]);
        c.LevelList = linspace(m,M,20);
        axis([0 L 0 H])
        pbaspect([L H 1])
    end
    
    pause();
    t = t + dt;
end