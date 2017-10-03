function [interpsol] = interp_2solutions_on_fespace(fespace1,sol1,fespace2,sol2,target_fespace)

nodes = target_fespace.nodes;

interpsol = zeros(size(nodes,1),1);

for i = 1:length(nodes)
    node = nodes(i,1:2);
    
    [I1,code1] = interpolate_in_point(fespace1,sol1,node(1),node(2));
    [I2,code2] = interpolate_in_point(fespace2,sol2,node(1),node(2));
    
    if (code1 == 0 && code2 == 1)
        interpsol(i) = I1;
    elseif (code2 == 0 && code1 == 1)
        interpsol(i) = I2;
    elseif (code1 == 0 && code2 == 0)
        interpsol(i) = I1;
    else
        error(['Point (',num2str(node(1)),',',num2str(node(2)),') is outside the domain']);
    end
end
