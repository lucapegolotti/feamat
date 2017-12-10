function export_vtk_fluid(sol,fespace_u,fespace_p,filename)
% Export velocity-pressure in vtk format by using vtkwrite
% input=
%           sol: vector of degrees of freedom of the solution
%           fespace: finite element space
%           filename: output filename
%

if (strcmp(fespace_u.mesh.type,'structured'))
    x = fespace_u.mesh.X(:,1);
    y = fespace_u.mesh.Y(1,:);
    [Y,X,Z] = meshgrid(x,y,0);
    % separate velocity and pressure
    n_nodes_u = size(fespace_u.nodes,1);
    n_nodes_p = size(fespace_p.nodes,1);
    
    u1 = sol(1:n_nodes_u);
    u2 = sol(n_nodes_u+1:2*n_nodes_u);
    p = sol(2*n_nodes_u+1:end);
    u1 = reshape(u1(1:size(fespace_u.mesh.vertices,1)),size(X,1),size(X,2));
    u2 = reshape(u2(1:size(fespace_u.mesh.vertices,1)),size(X,1),size(X,2));
    u3 = 0*u1;
    p = reshape(p(1:size(fespace_p.mesh.vertices,1)),size(X,1),size(X,2));
    
    vtkwrite(filename,'structured_grid',X,Y,Z,'scalars','pressure',p, ...
        'vectors','velocity', u1,u2,u3,'binary');
else
    error('Mesh type not supported!');
end