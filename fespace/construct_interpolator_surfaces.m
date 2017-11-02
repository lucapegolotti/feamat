function [R] = construct_interpolator_surfaces(fespace1,fespace2,index1,index2)
% construct the Lagrange interpolator from the boundary index1 of fespace1
% to the boundary index2 of fespace2 (needed for internodes approach)

warning('This function is only intended for use on structured meshes');

list1 = fespace1.boundary_nodes{index1};
list2 = fespace2.boundary_nodes{index2};

degree = fespace1.degree;

n1 = size(list2,1);
n2 = size(list1,1);

% we trasform the lists into vectors (one coordinate is always constant)
if (index1 == 1 || index1 == 3)
    list1 = list1(:,1);
else
    list1 = list1(:,2);
end

if (index2 == 1 || index2 == 3)
    list2 = list2(:,1);
else
    list2 = list2(:,2);
end

R = zeros(n1,n2);

l1 = list1(2)-list1(1);

iscontained = @(x1,x2,xp) (x1 <= xp)*(x2 >= xp);

% note: here we are using the fact that the mesh is structured (in
% particular the boundary are assumed to be traversed in the same
% direction)
if (strcmp(degree,'P1'))
    count = 1;
    x1 = list1(1);
    x2 = list1(2);
    i = 0;
    while (i <= n1)
        i = i+1;
        if (iscontained(x1,x2,list2(i)))
            l2 = list2(i)-x1;
            R(i,count) = 1-l2/l1;
            l2 = x2-list2(i);
            R(i,count+1) = 1-l2/l1;
        else
            count = count + 1;
            x1 = list1(count);
            x2 = list1(count+1);
            i = i-1;
        end
    end
elseif (strcmp(degree,'P2'))
    count = 1;
    
    % rearrange nodes such that middle points are between corner points
    R1 = zeros(n2);
    R2 = zeros(n1);
    
    corner = 1;
    middle = (n2-1)/2+2;
    row = 1;
    R1(row,corner) = 1;
    corner = corner + 1;
    while (row < n2)
        R1(row+1,middle) = 1;
        R1(row+2,corner) = 1;
        corner = corner + 1;
        row = row + 2;
        middle = middle + 1;
    end
    
    corner = 1;
    middle = (n1-1)/2+2;
    row = 1;
    R2(row,corner) = 1;
    corner = corner + 1;
    while (row < n1)
        R2(row+1,middle) = 1;
        R2(row+2,corner) = 1;
        corner = corner + 1;
        row = row + 2;
        middle = middle + 1;
    end
    
    list1 = R1*list1;
    list2 = R2*list2;
    
    x1 = list1(1);
    x2 = list1(2);
    x3 = list1(3);
    mat = [x1^2 x1 1;
        x2^2 x2 1;
        x3^2 x3 1];
    i = 0;
    while (i <= n1)
        i = i+1;
        if (iscontained(x1,x3,list2(i)))
            coeff = mat\[1;0;0];
            R(i,count) = coeff'*[list2(i)^2;list2(i);1];
            
            coeff = mat\[0;1;0];
            R(i,count+1) = coeff'*[list2(i)^2;list2(i);1];
            
            coeff = mat\[0;0;1];
            R(i,count+2) = coeff'*[list2(i)^2;list2(i);1];
            
        else
            count = count + 2;
            x1 = list1(count);
            x2 = list1(count+1);
            x3 = list1(count+2);
            mat = [x1^2 x1 1;
                   x2^2 x2 1;
                   x3^3 x3 1];
            i = i-1;
        end
        R
    end
    
    list1 = R1'*list1;
    list2 = R2'*list2;
else
    error(['Interpolation is not implemented for polynomials of type ',degree]);
end



end

