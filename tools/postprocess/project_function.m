function nodal_values = project_function(fespace, fun)
% Project analytical function onto finite element space
% input= 
%           fespace: finite element space
%           fun: anonymous function to project
%
% output= 
%           nodal_values: values of the function at the nodes of the mesh
%           
%

nodes = fespace.nodes;
n = size(nodes,1);
nodal_values = zeros(n,1);

for i = 1:n
    nodal_values(i) = fun(nodes(i,1:2)');
end