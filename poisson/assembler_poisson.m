function [A,b] = assembler_poisson(fespace,fun,mu,dirichlet_functions,neumann_functions)
tic 
bc_flags = fespace.bc;

thereisneumann = 1;

if (length(find(bc_flags)) == 4)
    thereisneumann = 0;
end

A = assemble_stiffness(mu,fespace);
b = assemble_rhs(fespace,fun);

if (thereisneumann)
   b = apply_neumann_bc(fespace,b,neumann_functions); 
end

% Apply Dirichlet boundary conditions
[A,b] = apply_dirichlet_bc(A,b,fespace,dirichlet_functions);

elapsed = toc;
disp(['Assembly of the system took = ', num2str(elapsed),' s']);
disp('------------------------------');


