function export_vtk_fluid(sol,filename)
% Export velocity-pressure in vtk format by using vtkwrite
% input=
%           sol: vector of degrees of freedom of the solution
%           fespace: finite element space
%           filename: output filename
%

fespace_u = sol.fespace_u;
fespace_p = sol.fespace_p;

if (strcmp(fespace_u.mesh.type,'structured'))
    x = fespace_u.mesh.X(:,1);
    y = fespace_u.mesh.Y(1,:);
    [Y,X,Z] = meshgrid(x,y,0);

    u1 = reshape(sol.u1(1:size(fespace_u.mesh.vertices,1)),size(X,1),size(X,2));
    u2 = reshape(sol.u2(1:size(fespace_u.mesh.vertices,1)),size(X,1),size(X,2));
    u3 = 0*u1;
    p = reshape(sol.p(1:size(fespace_p.mesh.vertices,1)),size(X,1),size(X,2));
    
    vtkwrite(filename,'structured_grid',X,Y,Z,'scalars','pressure',p, ...
        'vectors','velocity', u1,u2,u3,'binary');
else
    error('Mesh type not supported!');
end