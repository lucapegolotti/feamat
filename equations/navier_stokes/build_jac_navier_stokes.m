function [J] = build_jac_navier_stokes(A,u,fespace_u)
% Build jacobian of the Navier Stokes equations
%
% input=
%           A: anonymous function of the matrix of the system
%           u: convective velocity
%           fespace_u: finite element space of the velocity
%
% output= 
%           J: jacobian

n_nodes_u = size(fespace_u.nodes,1);

indices_u1 = 1:n_nodes_u;
indices_u2 = n_nodes_u+1:n_nodes_u*2;

J11 = assemble_jac_convective_term(fespace_u,u,1,1); 
J12 = assemble_jac_convective_term(fespace_u,u,1,2);
J21 = assemble_jac_convective_term(fespace_u,u,2,1);
J22 = assemble_jac_convective_term(fespace_u,u,2,2);

J = A(u);
J(indices_u1,indices_u1) = apply_dirichlet_bc_matrix(J(indices_u1,indices_u1) + J11,fespace_u,1);
J(indices_u1,indices_u2) = apply_dirichlet_bc_matrix(J(indices_u1,indices_u2) + J12,fespace_u,0);
J(indices_u2,indices_u1) = apply_dirichlet_bc_matrix(J(indices_u2,indices_u1) + J21,fespace_u,0);
J(indices_u2,indices_u2) = apply_dirichlet_bc_matrix(J(indices_u2,indices_u2) + J22,fespace_u,1);

J = sparse(J);