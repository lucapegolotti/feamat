function [A,b] = apply_dirichlet_bc(A,b,fespace,dirichlet_functions)

n_nodes = size(A,1);
dim2 = size(A,2);

nodes = fespace.nodes;
bc_flags = fespace.bc;

for i = 1:n_nodes
    if (nodes(i,3)~=0)
        if (bc_flags(nodes(i,3)))
            vd = dirichlet_functions(nodes(i,:)');
            A(i,:) = zeros(1,dim2);
            A(i,i) = 1;
            b(i) = vd(nodes(i,3));
        elseif (nodes(i,4) ~= 0 && bc_flags(nodes(i,4)))
            vd = dirichlet_functions(nodes(i,:)');
            A(i,:) = zeros(1,dim2);
            A(i,i) = 1;
            b(i) = vd(nodes(i,4));
        end
    end
end