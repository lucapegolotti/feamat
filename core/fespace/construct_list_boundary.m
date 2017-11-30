function [list] = construct_list_boundary(fespace,boundary_index)
nodes = fespace.nodes;
n_nodes = size(nodes,1);

list = zeros(n_nodes,2);
count = 0;
for i = 1:n_nodes
    if (nodes(i,end) == boundary_index || nodes(i,end-1) == boundary_index)
        count = count+1;
        list(count,:) = nodes(i,1:2);
    end
end

list = list(1:count,:);

end