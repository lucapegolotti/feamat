% common variables

%% Test 1: check Dirichlet indexes on structured meshes

mesh = create_mesh(0,0,1,1,10,10);
fespace = create_fespace(mesh,'P2',[0 1 0 0]);
n_nodes = size(fespace.nodes,1);

b = ones(n_nodes,1);
b = apply_dirichlet_bc_rhs(b,fespace,@(x) [0;0;0;0]);

for i = 1:n_nodes
    if (b(i) == 0)
        assert(fespace.nodes(i,1) == 1);
    end
    if (fespace.nodes(i,1) == 1)
        assert((b(i)) == 0);
    end
end

fespace = create_fespace(mesh,'P2',[1 0 0 0]);

b = ones(n_nodes,1);
b = apply_dirichlet_bc_rhs(b,fespace,@(x) [0;0;0;0]);

for i = 1:n_nodes
    if (b(i) == 0)
        assert(fespace.nodes(i,2) == 0);
    end
    if (fespace.nodes(i,2) == 0)
        assert((b(i)) == 0);
    end
end

fespace = create_fespace(mesh,'P2',[0 0 1 0]);

b = ones(n_nodes,1);
b = apply_dirichlet_bc_rhs(b,fespace,@(x) [0;0;0;0]);

for i = 1:n_nodes
    if (b(i) == 0)
        assert(fespace.nodes(i,2) == 1);
    end
    if (fespace.nodes(i,2) == 1)
        assert((b(i)) == 0);
    end
end

fespace = create_fespace(mesh,'P2',[0 0 1 0]);

b = ones(n_nodes,1);
b = apply_dirichlet_bc_rhs(b,fespace,@(x) [0;0;0;0]);

for i = 1:n_nodes
    if (b(i) == 0)
        assert(fespace.nodes(i,2) == 1);
    end
    if (fespace.nodes(i,2) == 1)
        assert((b(i)) == 0);
    end
end

fespace = create_fespace(mesh,'P2',[0 0 0 1]);

b = ones(n_nodes,1);
b = apply_dirichlet_bc_rhs(b,fespace,@(x) [0;0;0;0]);

for i = 1:n_nodes
    if (b(i) == 0)
        assert(fespace.nodes(i,1) == 0);
    end
    if (fespace.nodes(i,1) == 0)
        assert((b(i)) == 0);
    end
end



