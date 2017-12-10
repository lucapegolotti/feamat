function H = add_convective_term(H,u,fespace_u)
% Add convective term to diagonal blocks of stokes matrix.
% input=
%           H: stokes matrix
% output=
%           H: navier stokes matrix

n_nodes_u = size(fespace_u.nodes,1);

indices_u1 = 1:n_nodes_u;
indices_u2 = n_nodes_u+1:n_nodes_u*2;

C = apply_dirichlet_bc_matrix(assemble_convective_term(fespace_u,u),fespace_u,0);

H(indices_u1,indices_u1) = H(indices_u1,indices_u1) + C;
H(indices_u2,indices_u2) = H(indices_u2,indices_u2) + C;

