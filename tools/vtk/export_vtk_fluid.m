function export_vtk_fluid(sol,filename,varargin)
% Export velocity-pressure in vtk format by using vtkwrite
% input=
%           sol: structure of the solution 
%           filename: output filename
%           (optional)
%           'U' or 'P', to export only velocity or pressure
%

fespace_u = sol.fespace_u;
fespace_p = sol.fespace_p;

if (strcmp(fespace_u.mesh.type,'structured'))
    x = fespace_u.mesh.X(:,1);
    y = fespace_u.mesh.Y(1,:);
    [X,Y,Z] = meshgrid(x,y,0);

    u1 = reshape(sol.u1(1:size(fespace_u.mesh.vertices,1)),size(X,1),size(X,2));
    u2 = reshape(sol.u2(1:size(fespace_u.mesh.vertices,1)),size(X,1),size(X,2));
    u3 = 0*u1;
    p = reshape(sol.p(1:size(fespace_p.mesh.vertices,1)),size(X,1),size(X,2));
    
    if (nargin < 3)
    vtkwrite([filename,'.vtk'],'structured_grid',X',Y',Z','scalars','pressure',p, ...
        'vectors','velocity', u1,u2,u3,'binary');
    elseif (strcmp(varargin{1},'U'))
        vtkwrite([filename,'_velocity.vtk'],'structured_grid',X',Y',Z','vectors','velocity', u1,u2,u3,'binary');
    elseif (strcmp(varargin{1},'P'))
        vtkwrite([filename,'_pressure.vtk'],'structured_grid',X',Y',Z','scalars','pressure',p,'binary');
    else
        error('Export option not supported!');
    end
else
    mesh = fespace_u.mesh;
    x = mesh.vertices(:,1);
    y = mesh.vertices(:,2);

    n_vertices = length(x);

    data_title = 'fluid plot';

    if (nargin < 3)
        data_struct.type = 'scalar';
        data_struct.name = 'pressure';
        data_struct.data = sol.p(1:n_vertices);

        data_struct(2).type = 'vector';
        data_struct(2).name = 'velocity';
        data_struct(2).data = [sol.u1(1:n_vertices), ...
                               sol.u2(1:n_vertices), ...
                               sol.u1(1:n_vertices)*0];

        vtk_write_triangular_grid_and_data([filename,'.vtk'],data_title,[x y 0*x],mesh.elements(:,1:3),data_struct,false);
    elseif (strcmp(varargin{1},'U'))
        data_struct.type = 'vector';
        data_struct.name = 'velocity';
        data_struct.data = [sol.u1(1:n_vertices), ...
                            sol.u2(1:n_vertices), ...
                            sol.u1(1:n_vertices)*0];

        vtk_write_triangular_grid_and_data([filename,'_velocity.vtk'],data_title,[x y 0*x],mesh.elements(:,1:3),data_struct,false);
    elseif (strcmp(varargin{1},'P'))
        data_struct.type = 'scalar';
        data_struct.name = 'pressure';
        data_struct.data = sol.p(1:n_vertices);

        vtk_write_triangular_grid_and_data([filename,'_pressure.vtk'],data_title,[x y 0*x],mesh.elements(:,1:3),data_struct,false);
    else
        error('Export option not supported!');
    end

end