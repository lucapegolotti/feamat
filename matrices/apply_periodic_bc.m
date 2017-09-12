function [A,b] = apply_periodic_bc(A,b,fespace)

n_nodes = size(A,1);
dim2 = size(A,2);

nodes = fespace.nodes;
bc_flags = fespace.bc;

for i = 1:n_nodes
    if (nodes(i,3)~=0)
        % 2 is the flag of periodic bc
        if (bc_flags(nodes(i,3)) == 2)
            if (nodes(i,3) == 1)
                A(i,:) = zeros(1,dim2);
                A(i,i) = 1;
                A(i,fespace.periodic_b1(i)) = -1;
                b(i) = 0;
            elseif (nodes(i,3) == 2)
                A(i,:) = zeros(1,dim2);
                A(i,i) = 1;
                A(i,fespace.periodic_b2(i)) = -1;
                b(i) = 0;
            elseif (nodes(i,3) == 3)
                A(i,:) = zeros(1,dim2);
                A(i,i) = 1;
                A(i,fespace.periodic_b3(i)) = -1;
                b(i) = 0;
            elseif (nodes(i,3) == 4)
                A(i,:) = zeros(1,dim2);
                A(i,i) = 1;
                A(i,fespace.periodic_b4(i)) = -1;
                b(i) = 0;
            end  
        end
    end
     if (nodes(i,4)~=0)
        % 2 is the flag of periodic bc
        if (bc_flags(nodes(i,4)) == 2)
            if (nodes(i,4) == 1)
                A(i,:) = zeros(1,dim2);
                A(i,i) = 1;
                A(i,fespace.periodic_b1(i)) = -1;
                b(i) = 0;
            elseif (nodes(i,4) == 2)
                A(i,:) = zeros(1,dim2);
                A(i,i) = 1;
                A(i,fespace.periodic_b2(i)) = -1;
                b(i) = 0;
            elseif (nodes(i,4) == 3)
                A(i,:) = zeros(1,dim2);
                A(i,i) = 1;
                A(i,fespace.periodic_b3(i)) = -1;
                b(i) = 0;
            elseif (nodes(i,4) == 4)
                A(i,:) = zeros(1,dim2);
                A(i,i) = 1;
                A(i,fespace.periodic_b4(i)) = -1;
                b(i) = 0;
            end  
        end
    end
end