function b = apply_dirichlet_bc_rhs(b,fespace,dirichlet_functions)

n_nodes = size(b,1);

nodes = fespace.nodes;
bc_flags = fespace.bc;

for i = 1:n_nodes
    if (nodes(i,3)~=0)
        if (bc_flags(nodes(i,3)))
            vd = dirichlet_functions(nodes(i,:)');
            b(i) = vd(nodes(i,3));
        elseif (nodes(i,4) ~= 0 && bc_flags(nodes(i,4)))
            vd = dirichlet_functions(nodes(i,:)');
            b(i) = vd(nodes(i,4));
        end
    end
end
