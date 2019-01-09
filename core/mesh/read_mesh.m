function [mesh] = read_mesh(mesh_file,boundary_indicators)
% Reads mesh from a .msh file (for example, generate by gmsh)
% input=
%           mesh_file: path to msh.file
%           (optional)
%           boundary_indicators: vector valued anonymous function @(x)
%                                with number of components = number of
%                                boundaries. The vector takes values 1 in
%                                the component corresponding to a boundary
%                                if the point x belongs to the
%                                corresponding boundary, zero otherwise.
%                                It is necessary to provide this argument,
%                                if the .msh does not contain any physical
%                                entity to assign the boundary conditions.
%
% output=
%           mesh: mesh data structure


fid = fopen(mesh_file);

if (fid == -1)
    error('Could not open mesh file!');
end

if (nargin == 2)
    boundaries = cell(1,size(boundary_indicators([0;0]),1));
end

res = 1;
physical_names = 0;

allocated_elements_vertex = 0;

while(res ~= -1)
    res = fgets(fid);

    if (compare_string(res,'$PhysicalNames'))
        physical_names = 1;
        max_index = 0;
        res = fgets(fid);
        res = fgets(fid);
        while (~compare_string(res,'$EndPhysicalNames'))
            values = strsplit(res(1:end-1));
            index = str2num(values{2});
            max_index = max(max_index,index);
            res = fgets(fid);
        end
        boundaries = cell(1,max_index);
    end

    % start reading nodes
    if (compare_string(res,'$Nodes') == 1)

        if (~physical_names && nargin == 1)
            error(['When no physical elements are specified in the .msh file', ...
                   ' it is necessary to provide boundary indicators to the', ...
                   ' read_mesh function!']);
        end
        res = fgets(fid);
        num_elements = str2num(res(1:end-1));
        
        if (~allocated_elements_vertex)
            elements_containing_vertex = cell(num_elements,1);
            allocated_elements_vertex = 1;
        end
        % allocate vertices
        vertices = zeros(num_elements,4);
        count = 0;
        while (1)
            count = count + 1;
            res = fgets(fid);
            if (~compare_string(res,'$EndNodes'))
                numbers = strsplit(res(1:end-1));
                vertices(count,1) = str2double(numbers{2});
                vertices(count,2) = str2double(numbers{3});
            else
                break;
            end
        end
    end
    
    % start reading elements
    if (compare_string(res,'$Elements') == 1)
        res = fgets(fid);
        num_elements = str2num(res(1:end-1));
        % allocate elements
        elements = zeros(num_elements,4);
        count = 0;
        hmax = 0;
        while (1)
            res = fgets(fid);
            if (~compare_string(res,'$EndElements'))
                numbers = strsplit(res(1:end-1));
                if (strcmp(numbers{2},'1'))
                    v = str2double(numbers{6});
                    % the vertex has not been given a boundary flag yet
                    if (vertices(v,3) == 0 && ~physical_names)
                        aux = boundary_indicators(vertices(v,1:2));
                        boundary_indices = find(aux ~= 0);
                        if (length(boundary_indices) == 1)
                            vertices(v,3) = boundary_indices;
                        elseif (length(boundary_indices) == 2)
                            vertices(v,3:4) = boundary_indices;
                        elseif (length(boundary_indices) > 2)
                            error(['Boundary vertices cannot belong to more than 2', ... 
                                   ' boundaries!'])
                        end
                        flags = boundary_indices;
                    elseif (physical_names)
                        % if physical entities are present, we use the
                        % same flag
                        flag = str2num(numbers{4});
                        if (vertices(v,3) == 0)
                            vertices(v,3) = flag;
                        elseif (vertices(v,4) == 0 && vertices(v,3) ~= flag)
                            vertices(v,4) = flag;
                        elseif (vertices(v,3) ~= flag && ...
                                vertices(v,4) ~= flag)
                            error(['Boundary vertices cannot belong to more than 2', ... 
                                   ' boundaries!'])
                        end
                    end
                    v1 = v;
                    v = str2double(numbers{7});
                    % the vertex has not been given a boundary flag yet
                    if (vertices(v,3) == 0 && ~physical_names)
                        aux = boundary_indicators(vertices(v,1:2));
                        boundary_indices = find(aux ~= 0);
                        if (length(boundary_indices) == 1)
                            vertices(v,3) = boundary_indices;
                        elseif (length(boundary_indices) == 2)
                            vertices(v,3:4) = boundary_indices;
                        elseif (length(boundary_indices) > 2)
                            error(['Boundary vertices cannot belong to more than 2', ... 
                                   ' boundaries!'])
                        end
                        flag = intersect(vertices(v1,3:4),vertices(v,3:4));
                        if (flag(1) == 0)
                            flag = flag(2);
                        end
                     elseif (physical_names)
                        % if physical entities are present, we use the
                        % same flag
                        flag = str2num(numbers{4});
                        if (vertices(v,3) == 0)
                            vertices(v,3) = flag;
                        elseif (vertices(v,4) == 0 && vertices(v,3) ~= flag)
                            vertices(v,4) = flag;
                        elseif (vertices(v,3) ~= flag && ...
                                vertices(v,4) ~= flag)
                            error(['Boundary vertices cannot belong to more than 2', ... 
                                   ' boundaries!'])
                        end
                    end 
                    v2 = v;
                    boundaries{flag} = [boundaries{flag};v1 v2];
                end
                % we process the triangles
                if (strcmp(numbers{2},'2'))
                    count = count + 1;
                    num1 = str2double(numbers{6});
                    num2 = str2double(numbers{7});
                    num3 = str2double(numbers{8});

                    elements(count,1) = num1;
                    elements(count,2) = num2;
                    elements(count,3) = num3;

                    elements_containing_vertex{num1} = [elements_containing_vertex{num1};count];
                    elements_containing_vertex{num2} = [elements_containing_vertex{num2};count];
                    elements_containing_vertex{num3} = [elements_containing_vertex{num3};count];

                    hmax = max(hmax, norm(vertices(elements(count,2),1:2) - vertices(elements(count,1),1:2)));
                    hmax = max(hmax, norm(vertices(elements(count,3),1:2) - vertices(elements(count,2),1:2)));
                    hmax = max(hmax, norm(vertices(elements(count,1),1:2) - vertices(elements(count,3),1:2)));
                    if (vertices(elements(count,1),3) ~= 0)
                        elements(count,4) = vertices(elements(count,1),3);
                    elseif (vertices(elements(count,2),3) ~= 0)
                        elements(count,4) = vertices(elements(count,2),3);
                    elseif (vertices(elements(count,3),3) ~= 0)
                        elements(count,4) = vertices(elements(count,3),3);
                    end
                end
            else
                break;
            end
        end
        elements = elements(1:count,:);
    end
end
fclose(fid);

mesh.vertices = vertices;
mesh.elements_containing_vertex = elements_containing_vertex;
mesh.elements = elements;
mesh.boundaries = boundaries;
mesh.xp = min(vertices(:,1));
mesh.yp = min(vertices(:,2));
mesh.L = max(vertices(:,1)) - mesh.xp;
mesh.H = max(vertices(:,2)) - mesh.yp;
mesh.h = hmax;
mesh.triang = triangulation(mesh.elements(:,1:3),mesh.vertices(:,1),mesh.vertices(:,2));
mesh.type = 'unstructured';
    
    function result = compare_string(read_str,target)
        result = strcmp(read_str(1:min(length(read_str),length(target))),target) == 1;
    end

end
