function [index,x1,x2,x3,index_neigb] = find_element_containing_point(mesh,x_p)

x1 = [0;0];
x2 = [0;0];
x3 = [0;0];

if (isfield(mesh,'triang'))
    index = pointLocation(mesh.triang,x_p(:)');
    % point lies outside of triangulation
    if (index ~= index)
        index = -1;
        return
    end
    x1 = mesh.vertices(mesh.elements(index,1),1:2)';
    x2 = mesh.vertices(mesh.elements(index,2),1:2)';
    x3 = mesh.vertices(mesh.elements(index,3),1:2)';
    return
end

x = x_p(1);
y = x_p(2);

index = -1;
index_neigb = -1;

if (x < mesh.xp || x > mesh.xp+mesh.L)
    index = -1;
    return;
end
if (y < mesh.yp || y > mesh.yp+mesh.H)
    index = -1;
    return;
end

if (strcmp(mesh.type,'structured'))
    % find element that contains (x,y)
    xx = mesh.X(1:end-1,1);
    yy = mesh.Y(1,1:end-1);
    
    L = length(xx);
    
    indx = max(find(xx <= x));
    indy = max(find(yy <= y));
    
    elindex = 2*L*(indy-1) + indx;
    
    % check if (x,y) is in this element or the one above
    indv = mesh.elements(elindex,:);
    
    x1 = mesh.vertices(indv(1),1:2)';
    x2 = mesh.vertices(indv(2),1:2)';
    x3 = mesh.vertices(indv(3),1:2)';
    
    m1 = (x3(2)-x1(2))/(x3(1)-x1(1));
    m2 = (y-x1(2))/(x-x1(1));
    
    % then we chose the wrong element: choosing the one above
    if (m1 < m2)
        indv = mesh.elements(elindex+L,:);
        x1 = mesh.vertices(indv(1),1:2)';
        x2 = mesh.vertices(indv(2),1:2)';
        x3 = mesh.vertices(indv(3),1:2)';
        index = elindex+L;
        index_neigb = elindex;
        return;
    end
    index = elindex;
    index_neigb = elindex+L;
    return;
elseif (strcmp(mesh.type,'unstructured'))
    aux = mesh.vertices(:,1:2) - x_p(:)';
    distances = sqrt(aux(:,1).^2 + aux(:,2).^2);
    
    % search closest 2 vertices to input point
    [~,index1] = min(distances);
    distances(index1) = Inf;
    
    [~,index2] = min(distances);
    
    indices_elements = [mesh.elements_containing_vertex{index1}; ...
        mesh.elements_containing_vertex{index2}];
    
    found = 0;
    for i = 1:size(indices_elements,1)
        idx = indices_elements(i);
        x1 = mesh.vertices(mesh.elements(idx,1),1:2)';
        x2 = mesh.vertices(mesh.elements(idx,2),1:2)';
        x3 = mesh.vertices(mesh.elements(idx,3),1:2)';
        
        % check if the point is on a line
        if ((abs((norm(x1' - x_p) + norm(x2' - x_p))/norm(x1 - x2) - 1) < 1e-3) || ...
                (abs((norm(x2' - x_p) + norm(x3' - x_p))/norm(x2 - x3) - 1) < 1e-3) || ...
                (abs((norm(x3' - x_p) + norm(x1' - x_p))/norm(x3 - x1) - 1) < 1e-3))
            found = 1;
            break;
        else
            idx = indices_elements(i);
            x1 = mesh.vertices(mesh.elements(idx,1),1:2)';
            x2 = mesh.vertices(mesh.elements(idx,2),1:2)';
            x3 = mesh.vertices(mesh.elements(idx,3),1:2)';
            P12 = (x1-x2)'; P23 = (x2-x3)'; P31 = (x3-x1)';
            
            s1 = sign(det([P31;P23]))*sign(det([x3'-x_p(:)';P23]));
            s2 = sign(det([P12;P31]))*sign(det([x1'-x_p(:)';P31]));
            s3 = sign(det([P23;P12]))*sign(det([x2'-x_p(:)';P12]));
            is_inside =  s1 >= -1e-5 & s2 >= -1e-5 & s3 > -1e-5;
            
            if (is_inside)
                found = 1;
                break;
            end
        end
    end
        
    end
    if (found == 0)
        index = -1;
        return;
    end
    index = idx;
    return;
end