function [interpsol] = interp_on_fespace(fespace1,sol1,target_fespace)

nodes = target_fespace.nodes;

interpsol = zeros(size(nodes,1),1);

for i = 1:length(nodes)
    node = nodes(i,1:2);
    
    [I1,code1] = evaluate_fe_function(sol1,fespace1,node);
    
    if (code1 == 0)
        interpsol(i) = I1;
    else
        interpsol(i) = 0;
        %error(['Point (',num2str(node(1)),',',num2str(node(2)),') is outside the domain']);
    end
end
