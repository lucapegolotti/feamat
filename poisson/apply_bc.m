function [A,b] = apply_bc(A,b,nodes,bc_flags,dirichlet_functions)

n_nodes = size(A,1);

disp('Applying boundary conditions');
for i = 1:n_nodes
    if (nodes(i,3)~=0)
        if (bc_flags(nodes(i,3)))
            vd = dirichlet_functions(nodes(i,1),nodes(i,2));
            A(i,:) = zeros(1,n_nodes);
            A(i,i) = 1;
            b(i) = vd(nodes(i,3));
        end
    end
    
end
