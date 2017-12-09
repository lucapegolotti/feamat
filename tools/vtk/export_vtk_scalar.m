function export_vtk_scalar(sol,fespace,filename)
% Export solution in vtk format by using vtkwrite
% input=
%           sol: vector of degrees of freedom of the solution
%           fespace: finite element space
%           filename: output filename
%

if (strcmp(fespace.mesh.type,'structured'))
    x = fespace.mesh.X(:,1);
    y = fespace.mesh.Y(1,:);
    [Y,X,Z] = meshgrid(x,y,0);
    solgrid = reshape(sol(1:size(fespace.mesh.vertices,1)),size(X,1),size(X,2));
    vtkwrite(filename,'structured_grid',X,Y,Z,'scalars','u', solgrid,'binary')
else
    error('Mesh type not supported!');
end