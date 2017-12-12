function [x,u] = get_values_over_line(fespace,sol,samplepoints,coordinate,label)
xmin = fespace.mesh.xp;
ymin = fespace.mesh.yp;
xmax = fespace.mesh.xp + fespace.mesh.L;
ymax = fespace.mesh.yp + fespace.mesh.H;

u = zeros(samplepoints,1);
% parallel to the x axis
if (strcmp(label,'Xpar'))
    x = linspace(xmin,xmax,samplepoints);
    for j = 1:samplepoints
        u(j) = evaluate_fe_function(sol,fespace,[x(j);coordinate]);
    end
elseif (strcmp(label,'Ypar'))
    x = linspace(ymin,ymax,samplepoints);
    for j = 1:samplepoints
        u(j) = evaluate_fe_function(sol,fespace,[coordinate;x(j)]);
    end
end

end

