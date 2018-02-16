function streaml = plot_streamlines(sol,xstart,ystart,step)
% Plot streamlines of a fluid function
% input=
%           sol: structure containing the solution
%           xstart: vector of the x-coordinates from which streamlines should
%               be originated
%           ystart: vector of the y-coordinates from which streamlines should
%               be originated
%

if (~exist('step','var'))
    step = [];
end

fespace_u = sol.fespace_u;

if (strcmp(fespace_u.mesh.type,'structured'))
    n_vertices = size(fespace_u.mesh.vertices,1);
    
    x = fespace_u.mesh.X(:,1);
    y = fespace_u.mesh.Y(1,:);
    [X,Y] = meshgrid(x,y);
    
    n1 = size(fespace_u.mesh.X,1);
    n2 = size(fespace_u.mesh.X,2);
    U1 = reshape(sol.u1(1:n_vertices),n1,n2);
    U2 = reshape(sol.u2(1:n_vertices),n1,n2);
    
    hold on
    streams = mmstream2(X,Y,U1',U2',xstart,ystart,[],step);
    nstreams = size(streams,2);
    
    for i = 1:nstreams
        s = streams{i};
        s(s(:,1) == 0 & s(:,2) == 0,:) = [];
        streams{i} = s;
    end
    streaml = streamline(streams);
else
    error('Mesh type not supported!');
end