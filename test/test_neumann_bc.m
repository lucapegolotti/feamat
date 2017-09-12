clear allc
clc

mesh = create_mesh(1,1,2,2);

bc = [0 0 0 0];

fespace = create_fespace(mesh,'P1',bc);
 
neumann_bc =@(x) [1;0;1;0];

b = zeros(size(fespace.nodes,1),1);

b = apply_neumann_bc(fespace,b,neumann_bc)


