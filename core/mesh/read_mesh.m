function [mesh] = read_mesh(mesh_file,boundary_indicators)

fid = fopen(mesh_file);

if (fid == -1)
    error('Could not open mesh file!');
end

res = 1;
while(res ~= -1)
    res = fgets(fid);

    % start reading nodes
    if (strcmp(res(1:end-1),'$Nodes') == 1)
        res = fgets(fid);
        num_elements = str2num(res(1:end-1));

        % allocate vertices
        vertices = zeros(num_elements,4);
        count = 0;
        while (1)
            count = count + 1;
            res = fgets(fid);
            if (~strcmp(res(1:end-1),'$EndNodes'))
                numbers = strsplit(res(1:end-1));
                vertices(count,1) = str2double(numbers{2});
                vertices(count,2) = str2double(numbers{3});
            else
                break;
            end
        end
    end

    % start reading elements
    if (strcmp(res(1:end-1),'$Elements') == 1)
        res = fgets(fid);
        num_elements = str2num(res(1:end-1));

        % allocate elements
        elements = zeros(num_elements,4);
        count = 0;
        hmax = 0;
        while (1)
            res = fgets(fid);
            if (~strcmp(res(1:end-1),'$EndElements'))
                numbers = strsplit(res(1:end-1));

                if (strcmp(numbers{2},'1'))
                    v = str2double(numbers{6});
                    % the vertex has not been given a boundary flag yet
                    if (vertices(v,3) == 0)
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
                    end

                    v = str2double(numbers{7});
                    % the vertex has not been given a boundary flag yet
                    if (vertices(v,3) == 0)
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
                    end
                end
                % we process the triangles
                if (strcmp(numbers{2},'2'))
                    count = count + 1;
                    elements(count,1) = str2double(numbers{6});
                    elements(count,2) = str2double(numbers{7});
                    elements(count,3) = str2double(numbers{8});
                    
                    hmax = max(hmax, norm(vertices(elements(count,2),1:2) - vertices(elements(count,1),1:2)));
                    hmax = max(hmax, norm(vertices(elements(count,3),1:2) - vertices(elements(count,2),1:2)));
                    hmax = max(hmax, norm(vertices(elements(count,1),1:2) - vertices(elements(count,3),1:2)));
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
mesh.elements = elements;
mesh.xp = min(vertices(:,1));
mesh.yp = min(vertices(:,2));
mesh.L = max(vertices(:,1)) - mesh.xp;
mesh.H = max(vertices(:,1)) - mesh.yp;
mesh.h = hmax;
end

