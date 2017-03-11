function [X,Y,vertices, connectivity] = create_mesh(L,H,n_elements1,n_elements2)
    n_elements = n_elements1 * n_elements2 * 2;
    x = linspace(0,L,n_elements1+1);
    y = linspace(0,H,n_elements2+1);
    
    [X,Y] = meshgrid(x,y);
    X = X';
    Y = Y';
    vertices = [X(:) Y(:)];
    
    % last column represents group: 1: bottom, 2: right, 3: top, 4: left,
    %                               5: internal
    connectivity = zeros(n_elements,4);
    
    count = 0;
    for i = 1:n_elements2
       for j = 1:n_elements1
          count = count + 1;
          connectivity(count,:) = [j+(i-1)*(n_elements1+1) j+1+(i-1)*(n_elements1+1) j+1+i*(n_elements1+1) 5];
       end
       
       for j = 1:n_elements1
          count = count + 1;
          connectivity(count,:) = [j+(i-1)*(n_elements1+1) j+i*(n_elements1+1) j+i*(n_elements1+1)+1 5];
       end
    end
    
end