function [mesh] = create_mesh(xp,yp,L,H,n_elements1,n_elements2)
% Generates rectangular mesh.
% input= 
%        xp: x coordinate of the bottom left corner
%        yp: y coordinate of the bottom left corner
%        L: length
%        H: height
%        n_elements1: number of elements in the x direction
%        n_elements2: number of elements in the y direction
%
% output= X: x coordinates of the grid
%         Y: y coordinates of the grid
%         vertices: vertices of the mesh
%         connectivity: connectivity matrix
    
    n_elements = n_elements1 * n_elements2 * 2;
    x = linspace(xp,xp+L,n_elements1+1);
    y = linspace(yp,yp+H,n_elements2+1);
    
    [X,Y] = meshgrid(x,y);
    X = X';
    Y = Y';
    vertices = [X(:) Y(:) zeros(length(X(:)),1) zeros(length(X(:)),1)];
  
    % last column represents group: 1: bottom, 2: right, 3: top, 4: left,
    %                               0: internal
    elements = zeros(n_elements,4);
    
    disp(['Creating rectangular mesh with ', num2str(n_elements), ...
           ' elements (L = ', num2str(L), ', H = ', num2str(H),')']);
       
    tic
    
    count = 0;
    i = 1;
    for j = 1:n_elements1
        count = count + 1;
        if (j == n_elements1)
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+1+(i-1)*(n_elements1+1) j+1+i*(n_elements1+1) 2];
            vertices(j+(i-1)*(n_elements1+1),3) = 1;
            vertices(j+1+(i-1)*(n_elements1+1),3) = 1;
            vertices(j+1+(i-1)*(n_elements1+1),4) = 2;
            vertices(j+1+i*(n_elements1+1),3) = 2;
        else
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+1+(i-1)*(n_elements1+1) j+1+i*(n_elements1+1) 1];
            vertices(j+(i-1)*(n_elements1+1),3) = 1;
            vertices(j+1+(i-1)*(n_elements1+1),3) = 1;
        end
    end
       
    for j = 1:n_elements1
        count = count + 1;
        if (j == 1)
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+i*(n_elements1+1) j+i*(n_elements1+1)+1 4];
            vertices(j+(i-1)*(n_elements1+1),3) = 4;
            vertices(j+(i-1)*(n_elements1+1),4) = 1;
            vertices(j+i*(n_elements1+1),3) = 4;
        else
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+i*(n_elements1+1) j+i*(n_elements1+1)+1 0];
        end
    end
    
    for i = 2:n_elements2-1
       for j = 1:n_elements1
          count = count + 1;
          if (j == n_elements1)
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+1+(i-1)*(n_elements1+1) j+1+i*(n_elements1+1) 2];
            vertices(j+1+(i-1)*(n_elements1+1),3) = 2;
            vertices(j+1+i*(n_elements1+1),3) = 2;
          else
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+1+(i-1)*(n_elements1+1) j+1+i*(n_elements1+1) 0];
          end
       end
       
       for j = 1:n_elements1
          count = count + 1;
          if (j == 1)
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+i*(n_elements1+1) j+i*(n_elements1+1)+1 4];
            vertices(j+(i-1)*(n_elements1+1),3) = 4;
            vertices(j+i*(n_elements1+1),3) = 4;
          else
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+i*(n_elements1+1) j+i*(n_elements1+1)+1 0];
          end
       end
    end
    i = n_elements2;
    for j = 1:n_elements1
        count = count + 1;
        if (j == n_elements1)
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+1+(i-1)*(n_elements1+1) j+1+i*(n_elements1+1) 2];
            vertices(j+1+(i-1)*(n_elements1+1),3) = 2;
            vertices(j+1+i*(n_elements1+1),3) = 2;
            vertices(j+1+i*(n_elements1+1),4) = 2;
        else
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+1+(i-1)*(n_elements1+1) j+1+i*(n_elements1+1) 0];
        end
    end
       
    for j = 1:n_elements1
        count = count + 1;
        if (j == 1)
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+i*(n_elements1+1) j+i*(n_elements1+1)+1 4];
            vertices(j+(i-1)*(n_elements1+1),3) = 4;
            vertices(j+i*(n_elements1+1),3) = 4;
            vertices(j+i*(n_elements1+1),4) = 3;
        else
            elements(count,:) = [j+(i-1)*(n_elements1+1) j+i*(n_elements1+1) j+i*(n_elements1+1)+1 3];
            vertices(j+i*(n_elements1+1),3) = 3;
            vertices(j+i*(n_elements1+1)+1,3) = 3;
        end
    end
    
    mesh.vertices = vertices;
    mesh.elements = elements;
    mesh.xp = xp;
    mesh.yp = yp;
    mesh.X = X;
    mesh.Y = Y;
    mesh.L = L;
    mesh.H = H;
    
    elapsed = toc;
    disp(['Elapsed time = ', num2str(elapsed),' s']);
    disp('------------------------------');
    
end

