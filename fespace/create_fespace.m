function [fespace] = create_fespace(mesh,polydegree,bc_flags)

fespace.degree = 'P1';
% bcs are not enfornced at the level of the finite element space;
% they are enforced on the assembled matrix of the system
fespace.bc = bc_flags;

if (polydegree == 'P1')
    fespace.nodes = mesh.vertices;
    fespace.connectivity = mesh.elements;
    fespace.mesh = mesh;
    fespace.n_functions_per_element = 3;
    % coefficients of s and t and value for (0,0)
    fespace.functions = @(x) [1-x(1)-x(2);x(1);x(2)];
    fespace.grads = @(x) [-1 -1; 1 0; 0 1]';
else 
    error([polydegree, ' is not a valid polynomial degree!']);
end


end