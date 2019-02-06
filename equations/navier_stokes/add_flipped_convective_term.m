function H = add_flipped_convective_term(H,u,fespace_u)
% Add convective term to diagonal blocks of stokes matrix.
% input=
%           H: stokes matrix
% output=
%           H: navier stokes matrix

n_nodes_u = size(fespace_u.nodes,1);

indices_u1 = 1:n_nodes_u;
indices_u2 = n_nodes_u+1:n_nodes_u*2;

C_11 = apply_dirichlet_bc_matrix( assemble_flipped_convective_term(fespace_u,u,1,1), fespace_u, 0);
C_12 = apply_dirichlet_bc_matrix( assemble_flipped_convective_term(fespace_u,u,1,2), fespace_u, 0);
C_21 = apply_dirichlet_bc_matrix( assemble_flipped_convective_term(fespace_u,u,2,1), fespace_u, 0);
C_22 = apply_dirichlet_bc_matrix( assemble_flipped_convective_term(fespace_u,u,2,2), fespace_u, 0);

H(indices_u1,indices_u1) = H(indices_u1,indices_u1) + C_11;
H(indices_u1,indices_u2) = H(indices_u1,indices_u2) + C_12;
H(indices_u2,indices_u1) = H(indices_u2,indices_u1) + C_21;
H(indices_u2,indices_u2) = H(indices_u2,indices_u2) + C_22;

