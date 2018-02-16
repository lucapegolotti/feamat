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
    error('Mesh type not supported!');
end