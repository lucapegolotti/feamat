function [h]= plot_fe_function(vec,fespace,varargin)
% Plot finite element function as a surface (note: only value at vertices 
% of triangles are used for the visualization)
%
% input=
%           vec: vector of degrees of freedom
%           fespace: finite element space
%           (optional)
%           'contour' or 'contourf' to change visualization mode
%

n_vertices = size(fespace.mesh.vertices,1);

if (nargin == 3)
    error('Contour plots need the specification of the number of contourlines as 4th argument!');
end

if (strcmp(fespace.mesh.type,'structured'))
    n1 = size(fespace.mesh.X,1);
    n2 = size(fespace.mesh.X,2);

    if (nargin == 2)
        [h] = surf(fespace.mesh.X,fespace.mesh.Y,reshape(vec(1:n_vertices),n1,n2),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    else
        if (strcmp(varargin{1},'contourf'))
            [h] = contourf(fespace.mesh.X,fespace.mesh.Y,reshape(vec(1:n_vertices),n1,n2),varargin{2});
        elseif (strcmp(varargin{1},'contour'))
            [c,h] = contour(fespace.mesh.X,fespace.mesh.Y,reshape(vec(1:n_vertices),n1,n2),varargin{2},'ShowText','on');
            clabel(c,h,'labelspacing', 1000,'Fontsize',12);
        else
            error('Visualization option is not supported');
        end
    end
else
    if (nargin == 2)
        [h] = trisurf(fespace.mesh.elements(:,1:3),fespace.mesh.vertices(:,1), ...
                      fespace.mesh.vertices(:,2),vec(1:n_vertices));
    else
        if (strcmp(varargin{1},'contourf'))
            [h] = tricontf(fespace.mesh.vertices(:,1),fespace.mesh.vertices(:,2), ...
                           fespace.mesh.elements(:,1:3),vec(1:n_vertices),varargin{2});
        elseif (strcmp(varargin{1},'contour'))
            error(['Simple contour plot is not implemented for plots on', ...
                   ' unstructured meshes!']);
        else
            error('Visualization option is not supported');
        end
    end
end
