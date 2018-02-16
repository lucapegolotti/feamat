function plot_closed_streamlines(sol,xstart,ystart,xcenter,ycenter,color,res,step)
% Plot streamlines of a fluid function
% input=
%           sol: structure containing the solution
%           xstart: vector of the x-coordinates from which streamlines should
%               be originated
%           ystart: vector of the y-coordinates from which streamlines should
%               be originated
%           xcenter: x-coordinate of the center (any point inside the loop)
%           ycenter: y-coordinate of the center (any point inside the loop)
%           color: color of the lines
%           res: resolution of the lines, i.e. number of points between
%               each plotted point
%           (optional)
%           step: timestep for the integration
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
        angle = atan((s(:,2)-ycenter)./(s(:,1)-xcenter));
        first = 0;
        
        nangl = length(angle);
        count = 0;
        if (angle(2) - angle(1) > 0)
            for j = 1:nangl-1
                if (angle(j+1) > angle(1) && angle(j) <= angle(1))
                    count = count + 1;
                    if (count == 3)
                        first = j;
                        break;
                    end
                end
            end
        else
            for j = 1:nangl-1
                if (angle(j+1) < angle(1) && angle(j) >= angle(1))
                    count = count + 1;
                    if (count == 3)
                        first = j;
                        break;
                    end
                end
            end
        end
        if (first ~= 0)
            s = s(1:first,:);
        end
        s = [s(1:res:end,:);s(end,:)];
        plot(s(:,1),s(:,2),'Color',color,'Linewidth',1)
    end
else
    error('Mesh type not supported!');
end