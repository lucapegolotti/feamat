function [A,rhs] = assembler_diffusion_transport_advection(fespace,fun,mu,b,c,dirichlet_functions,neumann_functions)
% builds matrices and rhs for problem of the form: 
% -mu div(grad(u)) +  b grad u + c u = fun

tic 
bc_flags = fespace.bc;

thereisneumann = 1;

if (length(find(bc_flags)) == 4)
    thereisneumann = 0;
end

A = assemble_stiffness(mu,fespace) + assemble_transport(b,fespace) + ...
    assemble_advection(c,fespace);
rhs = assemble_rhs(fespace,fun);

if (thereisneumann)
   rhs = apply_neumann_bc(fespace,b,neumann_functions); 
end

% Apply Dirichlet boundary conditions
[A,rhs] = apply_dirichlet_bc(A,rhs,fespace,dirichlet_functions);

elapsed = toc;
disp(['Assembly of the system took = ', num2str(elapsed),' s']);
disp('------------------------------');


