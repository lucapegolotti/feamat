function b = apply_neumann_bc(b,fespace,neumann_functions,varargin)
% Apply Neumann boundary conditions to matrix to rhs
%
% input=
%           b: right handside
%           fespace: finite element space
%           neumann_functions: Neumann data
%           (optional)
%           order of gauss points (default = 2)
%
% output=
%           b: modified right handside
n_gauss1d = 2;
if (nargin == 4)
    n_gauss1d = varargin{1};
end

connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nlocalfunctions = fespace.n_functions_per_element;
bc_flags = fespace.bc;

n_elements = size(connectivity,1);

[gp1d,weights1d,~] = gauss_points1D(n_gauss1d);

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
            if (common_boundary ~= 0 && bc_flags(common_boundary)==0)
                transfgp_1d_to_2d = @(x) [1;0]/2*x + [1;0]/2;
                transf1d = @(x) (x2_edge - x1_edge)/2*x + (x1_edge + x2_edge)/2;
                dtransf = norm((x2_edge - x1_edge)/2);
                
                for j = 1:n_gauss1d
                    functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
                    nf = neumann_functions(transf1d(gp1d(j)));
                    for k = 1:nlocalfunctions
                        b(indices(k)) = b(indices(k)) + dtransf*nf(common_boundary)*functions(k)*weights1d(j);
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
            if (common_boundary ~= 0 && bc_flags(common_boundary)==0)
                transfgp_1d_to_2d = @(x) [-1;1]/2*x + [1;1]/2;
                transf1d = @(x) (x2_edge - x1_edge)/2*x + (x1_edge + x2_edge)/2;
                dtransf = norm((x2_edge - x1_edge)/2);
                
                for j = 1:n_gauss1d
                    functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
                    nf = neumann_functions(transf1d(gp1d(j)));
                    for k = 1:nlocalfunctions
                        b(indices(k)) = b(indices(k)) + dtransf*nf(common_boundary)*functions(k)*weights1d(j);
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
            if (common_boundary ~= 0 && bc_flags(common_boundary)==0)
                transfgp_1d_to_2d = @(x) [0;-1]/2*x + [0;1]/2;
                transf1d = @(x) (x2_edge - x1_edge)/2*x + (x1_edge + x2_edge)/2;
                dtransf = norm((x2_edge - x1_edge)/2);
                
                for j = 1:n_gauss1d
                    functions = fespace.functions(transfgp_1d_to_2d(gp1d(j)));
                    nf = neumann_functions(transf1d(gp1d(j)));
                    for k = 1:nlocalfunctions
                        b(indices(k)) = b(indices(k)) + dtransf*nf(common_boundary)*functions(k)*weights1d(j);
                    end
                end
            end
        end
    end
end
