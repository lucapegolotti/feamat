function plot_solution_vp(fespace_u,fespace_p,sol, what)

if (what == 'U')
    n_nodes_u = size(fespace_u.nodes,1);

    u1 = sol(1:n_nodes_u);
    u2 = sol(n_nodes_u+1:2*n_nodes_u);

    n_vertices = size(fespace_u.mesh.vertices,1);
    n1 = size(fespace_u.mesh.X,1);
    n2 = size(fespace_u.mesh.X,2);
    U1 = reshape(u1(1:n_vertices),n1,n2);
    U2 = reshape(u2(1:n_vertices),n1,n2);

    N = sqrt(U1.^2+U2.^2);
    m = min(min(N));
    M = max(max(N));


    [~,c] = contourf(fespace_u.mesh.X,fespace_u.mesh.Y,N);
    c.LineStyle = 'none';
    c.LevelStep = (M-m)/20;
    shading interp
    colorbar
    hold on

    q = quiver(fespace_u.mesh.X,fespace_u.mesh.Y,U1./N,U2./N);
    q.Color = 'black';
    q.LineWidth = 1;
    q.AutoScaleFactor = 0.5;
    hold off
elseif (what == 'P')
    n_nodes_u = size(fespace_u.nodes,1);

    p = sol(n_nodes_u*2+1:end);

    n_vertices = size(fespace_u.mesh.vertices,1);
    n1 = size(fespace_u.mesh.X,1);
    n2 = size(fespace_u.mesh.X,2);
   
   
    m = min(p);
    M = max(p);


    [~,c] = contourf(fespace_u.mesh.X,fespace_u.mesh.Y,reshape(p,n1,n2));
    c.LineStyle = 'none';
    c.LevelStep = abs(M-m)/20;
    colorbar    
end