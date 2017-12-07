function [fespace] = create_fespace(mesh,polydegree,bc_flags)
% Generates finite element space
% input= 
%           mesh: mesh data structure (created with create_mesh)
%           polydegree: string for polynomial degree. Currently supported: 
%                       'P1' and 'P2'
%           bc_flags: vector of flag for boundary conditions. We expect 4
%                     flags, ordered from the bottom edge in
%                     counterclockwise order. If a flag is equal to 0, then
%                     we consider Neumann boundary conditions. If a flag is
%                     equal to 1, we consider Dirichlet boundary conditions
% output= 
%           fespace: a fespace datastructure

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
    % compute nodes inside the triangles
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

end