% common variables

tol = 1e-12;

mesh = create_mesh(0,0,1,1,1,1);

fespaceP1 = create_fespace(mesh,'P1',[0 0 0 0]);
fespaceP2 = create_fespace(mesh,'P2',[0 0 0 0]);

%% Test 1: all fields assigned correctly

assert(strcmp(fespaceP1.degree,'P1'))
assert(strcmp(fespaceP2.degree,'P2'))

assert(norm(fespaceP1.nodes-mesh.vertices) < tol)
assert(norm(fespaceP1.connectivity-mesh.elements) < tol)
assert(fespaceP1.n_functions_per_element == 3)

expected_nodesP2 = [0 0 1 4;
                    1 0 1 2;
                    0 1 4 3;
                    1 1 2 3;
                    0.5 0 1 0;
                    1 0.5 2 0;
                    0.5 0.5 0 0;
                    0 0.5 4 0;
                    0.5 1.0 3 0];
                
assert(norm(fespaceP2.nodes-expected_nodesP2) < tol)

%% Test 2: value of basis functions is correct for P1
assert(norm(fespaceP1.functions([0;0])-[1;0;0]) < tol)
assert(norm(fespaceP1.functions([1;0])-[0;1;0]) < tol)
assert(norm(fespaceP1.functions([0;1])-[0;0;1]) < tol)

assert(norm(fespaceP2.functions([0;0])-[1;0;0;0;0;0]) < tol)
assert(norm(fespaceP2.functions([1;0])-[0;1;0;0;0;0]) < tol)
assert(norm(fespaceP2.functions([0;1])-[0;0;1;0;0;0]) < tol)
assert(norm(fespaceP2.functions([0.5;0])-[0;0;0;1;0;0]) < tol)
assert(norm(fespaceP2.functions([0.5;0.5])-[0;0;0;0;1;0]) < tol)
assert(norm(fespaceP2.functions([0;0.5])-[0;0;0;0;0;1]) < tol)



