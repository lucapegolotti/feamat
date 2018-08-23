function [M] = assemble_boundary_mass(fespace,flag)

n_gauss1d = 2;
if (nargin == 4)
    n_gauss1d = varargin{1};
end

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nlocalfunctions = fespace.n_functions_per_element;
bc_flags = fespace.bc;

n_elements = size(connectivity,1);
n_nodes = size(fespace.nodes,1);

[gp1d,weights1d,~] = gauss_points1D(n_gauss1d);

n_functions = fespace.n_functions_per_element;
n_functionsqr = n_functions^2;

elements_M = zeros(nlocalfunctions*n_elements,1);
indices_i = zeros(nlocalfunctions*n_elements,1);
indices_j = zeros(nlocalfunctions*n_elements,1);

for i = 1:n_elements
    indices = connectivity(i,1:end-1);
    x1 = vertices(indices(1),:)';
    x2 = vertices(indices(2),:)';
    x3 = vertices(indices(3),:)';
    
    % if boundary element, find vertices that are on edges
    if (connectivity(i,end) > 0)
        if (vertices(indices(1),3) > 0 && vertices(indices(2),3) > 0)
            x1_edge = x1;
            x2_edge = x2;
            common_boundary = 0;
            if (x1_edge(3) == x2_edge(3))
                common_boundary = x1_edge(3);
            elseif (x1_edge(4) == x2_edge(3))
                common_boundary = x1_edge(4);
            elseif (x1_edge(3) == x2_edge(4))
                common_boundary = x1_edge(3);
            elseif (x1_edge(4) == x2_edge(4))
                common_boundary = x1_edge(4);
            end
            
            x1_edge = x1_edge(1:2);
            x2_edge = x2_edge(1:2);
            
            
            % only continue if common boundary and common boundary == flag
            if (common_boundary ~= 0 && common_boundary == flag)
                [I1,I2] = meshgrid(indices,indices);
            
                currindices = (i-1)*n_functionsqr+1:n_functionsqr*i;
                indices_i(currindices) = I1(:);
                indices_j(currindices) = I2(:);

                new_elements = zeros(size(I1,1)^2,1);

                transfgp_1d_to_2d = @(x) [1;0]/2*x + [1;0]/2;
                dtransf = norm((x2_edge - x1_edge)/2);
                
                for j = 1:n_gauss1d
                    functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
                    matrix_elements = dtransf*(functions*functions')*weights1d(j);
                    new_elements = new_elements + matrix_elements(:);
                end
                elements_M(currindices) = new_elements;
            end
        end
        
        if (vertices(indices(2),3) > 0 && vertices(indices(3),3) > 0)
            x1_edge = x2;
            x2_edge = x3;
            common_boundary = 0;
            if (x1_edge(3) == x2_edge(3))
                common_boundary = x1_edge(3);
            elseif (x1_edge(4) == x2_edge(3))
                common_boundary = x1_edge(4);
            elseif (x1_edge(3) == x2_edge(4))
                common_boundary = x1_edge(3);
            elseif (x1_edge(4) == x2_edge(4))
                common_boundary = x1_edge(4);
            end
            
            x1_edge = x1_edge(1:2);
            x2_edge = x2_edge(1:2);         
            
            % only continue if common boundary and common boundary == flag
            if (common_boundary ~= 0 && common_boundary == flag)
                [I1,I2] = meshgrid(indices,indices);
            
                currindices = (i-1)*n_functionsqr+1:n_functionsqr*i;
                indices_i(currindices) = I1(:);
                indices_j(currindices) = I2(:);

                new_elements = zeros(size(I1,1)^2,1);
                transfgp_1d_to_2d = @(x) [-1;1]/2*x + [1;1]/2;
                dtransf = norm((x2_edge - x1_edge)/2);
                
                for j = 1:n_gauss1d
                    functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
                    matrix_elements = dtransf*(functions*functions')*weights1d(j);
                    new_elements = new_elements + matrix_elements(:);
                end
                elements_M(currindices) = new_elements;
            end
        end
        
        if (vertices(indices(3),3) > 0 && vertices(indices(1),3) > 0)
            x1_edge = x3;
            x2_edge = x1;
            common_boundary = 0;
            if (x1_edge(3) == x2_edge(3))
                common_boundary = x1_edge(3);
            elseif (x1_edge(4) == x2_edge(3))
                common_boundary = x1_edge(4);
            elseif (x1_edge(3) == x2_edge(4))
                common_boundary = x1_edge(3);
            elseif (x1_edge(4) == x2_edge(4))
                common_boundary = x1_edge(4);
            end
            
            x1_edge = x1_edge(1:2);
            x2_edge = x2_edge(1:2);
            
            % only continue if common boundary and common boundary == flag
            if (common_boundary ~= 0 && common_boundary == flag)          
                [I1,I2] = meshgrid(indices,indices);
            
                currindices = (i-1)*n_functionsqr+1:n_functionsqr*i;
                indices_i(currindices) = I1(:);
                indices_j(currindices) = I2(:);

                new_elements = zeros(size(I1,1)^2,1);
                transfgp_1d_to_2d = @(x) [0;-1]/2*x + [0;1]/2;
                dtransf = norm((x2_edge - x1_edge)/2);
                
                for j = 1:n_gauss1d
                    functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
                    matrix_elements = dtransf*(functions*functions')*weights1d(j);
                    new_elements = new_elements + matrix_elements(:);
                end
                elements_M(currindices) = new_elements;
            end
        end
    end
end

real_indices = find(indices_i ~= 0);

M = sparse(indices_i(real_indices),indices_j(real_indices),elements_M(real_indices),n_nodes,n_nodes);
