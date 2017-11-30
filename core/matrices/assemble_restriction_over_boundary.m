function [R] = assemble_restriction_over_boundary(fespace,boundary_index)

list = fespace.boundary_nodes{boundary_index};
nodes = fespace.nodes;

n_nodes = size(nodes,1);

R = zeros(size(list,1),n_nodes);

count = 0;
for i = 1:n_nodes
    if (nodes(i,end-1) == boundary_index || nodes(i,end) == boundary_index)
        count = count + 1
        R(count,i) = 1;
    end
    size(R)
end

%R = sparse(R);
end

