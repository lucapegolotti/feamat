function [fespace] = create_fespace(mesh,polydegree,bc_flags)

fespace.degree = polydegree;
% bcs are not enforced at the level of the finite element space;
% they are enforced on the assembled matrix of the system
fespace.bc = bc_flags;

if (polydegree == 'P1')
    fespace.nodes = mesh.vertices;
    fespace.connectivity = mesh.elements;
    fespace.mesh = mesh;
    fespace.n_functions_per_element = 3;
    fespace.functions = @(x) [1-x(1)-x(2);x(1);x(2)];
    fespace.grads = @(x) [-1 -1; 1 0; 0 1]';
elseif (polydegree == 'P2')
    n_elements = size(mesh.elements,1);
    n_vertices = size(mesh.vertices,1);
    
    aux = sparse(n_vertices,n_vertices);
    
    nodes = [mesh.vertices; zeros(size(mesh.vertices))];
    fespace.connectivity = zeros(n_elements,7);
    count = n_vertices + 1;
    for i = 1:n_elements
        
        % select indices and coordinates of the ith element
        indices = mesh.elements(i,:);
        X = [mesh.vertices(indices(1),:);mesh.vertices(indices(2),:);mesh.vertices(indices(3),:)];
        new_connectivity = zeros(1,6);
        for k = 1:3
            index1 = k;
            index2 = mod(k,3)+1;
            x1 = X(index1,:);
            x2 = X(index2,:);
            if (aux(indices(index1),indices(index2)) == 0)
                aux(indices(index1),indices(index2)) = count;
                bc = 0;
                for j = 3:4
                    if (x1(j) ~= 0 && (x1(j) == x2(3) || x1(j) == x2(4)))
                        bc = x1(j);
                    end
                end
                % note: this fails if the two nodes belong to different boundaries
                nodes(count,:) = [(x1(1:2) + x2(1:2))/2 bc 0];
                count = count + 1;
            end
            new_connectivity([k 3+k]) = [indices(index1) aux(indices(index1),indices(index2))];
        end
        fespace.connectivity(i,:) = [new_connectivity mesh.elements(i,4)];
    end
    fespace.nodes = nodes;
    fespace.mesh = mesh;
    fespace.n_functions_per_element = 6;
    c = [2 2 4 -3 -3 1; 2 0 0 -1 0 0; 0 2 0 0 -1 0; -4 0 -4 4 0 0; 0 0 4 0 0 0; 0 -4 -4 0 4 0];
    fespace.functions = @(x) c * [x(1)^2; x(2)^2; x(1)*x(2);x(1);x(2);1];
    fespace.grads = @(x) [4*x(1)+4*x(2)-3 4*x(1)+4*x(2)-3; ...
        4*x(1)-1 0; ...
        0 4*x(2)-1; ...
        -8*x(1)-4*x(2)+4 -4*x(1); ...
        4*x(2) 4*x(1); ...
        -4*x(2) -8*x(2)-4*x(1)+4 ]';
    
else
    error([polydegree, ' is not a valid polynomial degree!']);
end

indices_b1 = [];
indices_b2 = [];
indices_b3 = [];
indices_b4 = [];

nnodes = length(fespace.nodes);

for i = 1:nnodes
    b = fespace.nodes(i,3);
    if (b == 1)
        indices_b1 = [indices_b1; i];
    elseif (b == 2)
        indices_b2 = [indices_b2; i];
    elseif (b == 3)
        indices_b3 = [indices_b3; i];
    elseif (b == 4)
        indices_b4 = [indices_b4; i];
    end
    
    b = fespace.nodes(i,4);
    if (b == 1)
        indices_b1 = [indices_b1; i];
    elseif (b == 2)
        indices_b2 = [indices_b2; i];
    elseif (b == 3)
        indices_b3 = [indices_b3; i];
    elseif (b == 4)
        indices_b4 = [indices_b4; i];
    end
end
fespace.periodic_b1 = sparse(nnodes,1);
fespace.periodic_b1(indices_b1) = indices_b3;

fespace.periodic_b3 = sparse(nnodes,1);
fespace.periodic_b3(indices_b3) = indices_b1;

fespace.periodic_b2 = sparse(nnodes,1);
fespace.periodic_b2(indices_b2) = indices_b4;

fespace.periodic_b4 = sparse(nnodes,1);
fespace.periodic_b4(indices_b4) = indices_b2;

% check if correspondence is correct
for i = 1:length(indices_b1)
    if (fespace.nodes(indices_b1(i),1) ~= fespace.nodes(indices_b3(i),1))
        error('Error in periodic boundary detection!');
    end
end

for i = 1:length(indices_b2)
    if (fespace.nodes(indices_b2(i),2) ~= fespace.nodes(indices_b4(i),2))
        error('Error in periodic boundary detection!');
    end
end

fespace.boundary_nodes = {construct_list_boundary(fespace,1), ...
                          construct_list_boundary(fespace,2), ...
                          construct_list_boundary(fespace,3), ...
                          construct_list_boundary(fespace,4)};

end