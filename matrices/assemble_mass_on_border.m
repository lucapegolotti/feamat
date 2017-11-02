function [Mg] = assemble_mass_on_border(fespace,boundary_index)

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;
nlocalfunctions = fespace.n_functions_per_element;
bc_flags = fespace.bc;

n_elements = size(connectivity,1);

[gp1d,weights1d,n_gauss1d] = gauss_points1D(3);

Mg = sparse(size(nodes,1),size(nodes,1));

for i = 1:n_elements
    indices = connectivity(i,:);
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
            
            % only continue if common boundary and if neumann_conditions
            % must be applied (bc_flag == 0)
            if (common_boundary == boundary_index)
                transfgp_1d_to_2d = @(x) [1;0]/2*x + [1;0]/2;
                dtransf = norm((x2_edge - x1_edge)/2);
                
                for j = 1:n_gauss1d
                    functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
                    for k = 1:nlocalfunctions
                        elements = dtransf*(functions*functions')*weights1d(j);
                        Mg(indices(1:end-1),indices(1:end-1)) = Mg(indices(1:end-1),indices(1:end-1)) + elements;
                    end
                end
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
            
            % only continue if common boundary and if neumann_conditions
            % must be applied (bc_flag == 0)
            if (common_boundary == boundary_index)
                transfgp_1d_to_2d = @(x) [-1;1]/2*x + [1;1]/2;
                dtransf = norm((x2_edge - x1_edge)/2);
                
                for j = 1:n_gauss1d
                    functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
                    for k = 1:nlocalfunctions
                        elements = dtransf*(functions*functions')*weights1d(j);
                        Mg(indices(1:end-1),indices(1:end-1)) = Mg(indices(1:end-1),indices(1:end-1)) + elements;
                    end
                end
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
            
            % only continue if common boundary and if neumann_conditions
            % must be applied (bc_flag == 0)
            if (common_boundary == boundary_index)
                transfgp_1d_to_2d = @(x) [0;-1]/2*x + [0;1]/2;
                dtransf = norm((x2_edge - x1_edge)/2);
                
                for j = 1:n_gauss1d
                    functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
                    for k = 1:nlocalfunctions
                        elements = dtransf*(functions*functions')*weights1d(j);
                        Mg(indices(1:end-1),indices(1:end-1)) = Mg(indices(1:end-1),indices(1:end-1)) + elements;
                    end
                end
            end
        end
    end
end

Mg(~any(Mg,2),:) = [];
Mg(:,~any(Mg,1)) = [];


