function indices = find_dirichlet_indices(fespace)
% Find dirichlet indices
%
% input=
%           fespace: finite element space
%                                
% output=
%           indices: dirichlet indices
nodes = fespace.nodes;
bc_flags = fespace.bc;

n_nodes = size(nodes,1);

indices = [];
for i = 1:n_nodes
    if (nodes(i,3)~=0 && bc_flags(nodes(i,3)) || ...
        nodes(i,4) ~= 0 && bc_flags(nodes(i,4)))
        indices = [indices;i];
    end
end