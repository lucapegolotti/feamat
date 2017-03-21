function [A,b] = apply_bc(A,b,fespace,dirichlet_functions)

n_nodes = size(A,1);

nodes = fespace.nodes;
bc_flags = fespace.bc;

disp('Applying boundary conditions');
tic 
for i = 1:n_nodes
    if (nodes(i,3)~=0)
        if (bc_flags(nodes(i,3)))
            vd = dirichlet_functions(nodes(i,:)');
            A(i,:) = zeros(1,n_nodes);
            A(i,i) = 1;
            b(i) = vd(nodes(i,3));
        end
    end
end
    
elapsed = toc;
disp(['Elapsed time = ', num2str(elapsed),' s']);
disp('------------------------------');