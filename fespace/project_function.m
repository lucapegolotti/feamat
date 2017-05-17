function nodal_values = project_function(fespace, u)
nodes = fespace.nodes;
n = size(nodes,1);
nodal_values = zeros(n,1);

for i = 1:n
    nodal_values(i) = u(nodes(1:2)');
end