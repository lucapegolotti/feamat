function [v,vt] = assemble_vector_lagrange_multiplier_boundary(fespace,boundary_index)

v = zeros(size(fespace.nodes,1),1);
artificial_bc = [1 1 1 1];
artificial_bc(boundary_index) = 0;
truebc = fespace.bc;
fespace.bc = artificial_bc;

v = apply_neumann_bc(fespace,v,@(x) [1;1;1;1]); % attention, this may work only with scalar problems
vt = v';
v = apply_dirichlet_bc_rhs(v,fespace,@(x) [0;0;0;0]);