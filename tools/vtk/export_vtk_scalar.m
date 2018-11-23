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
    mesh = fespace.mesh;
    x = mesh.vertices(:,1);
    y = mesh.vertices(:,2);
    n_vertices = length(x);
    data_title = 'scalar plot';
    data_struct.type = 'scalar';
    data_struct.name = 'u';
    data_struct.data = sol;
    vtk_write_triangular_grid_and_data([filename,'.vtk'],data_title,[x y 0*x],mesh.elements(:,1:3),data_struct,false);
end