function [b] = assemble_rhs(fespace,fun)
% Assemble right handside
% input=
%           fespace: finite element space
%           fun: anonymous function or scalar (if fun == 0) of the right handside
%                If function, fun should accept matrices with two rows 
%                (first row x values, second row y values) to allow
%                vectorization when using structured meshes
%
% output=
%           b: right handside


connectivity = fespace.connectivity;
vertices = fespace.mesh.vertices;
nodes = fespace.nodes;
nlocalfunctions = fespace.n_functions_per_element;

n_elements = size(connectivity,1);
n_nodes = size(nodes,1);

if (~isa(fun,'function_handle'))
    if (fun ~= 0)
        error('Function must be a function handler or, if scalar, equal to zero!');
    end
    b = zeros(n_nodes,1);
    return;
end

n_gauss = 6;
[gp,weights,~] = gauss_points2D(n_gauss);

b = zeros(n_nodes,1);

if (~strcmp(fespace.mesh.type,'structured'))
    for i = 1:n_elements
        indices = connectivity(i,1:end-1);
        x1 = vertices(indices(1),1:2)';
        x2 = vertices(indices(2),1:2)';
        x3 = vertices(indices(3),1:2)';
        
        mattransf = [x2-x1 x3-x1];
        
        % transformation from parametric to physical
        transf = @(x) mattransf*x + x1;
        dettransf = abs(det(mattransf));
        
        for j = 1:n_gauss
            functions = fespace.functions(gp(:,j));
            for k = 1:nlocalfunctions
                b(indices(k)) = b(indices(k)) + dettransf*fun(transf(gp(:,j)))*functions(k)*weights(j)/2;
            end
        end
    end
else
    [fespace,~,~] = add_members_structured_meshes(fespace,n_gauss);
    
    rhs_elements1 = zeros(nlocalfunctions,n_gauss);
    rhs_elements2 = zeros(nlocalfunctions,n_gauss);
    for i = 1:n_gauss
        rhs_elements1(:,i) = fespace.dettransf1* ...
            fespace.transffuns1{i}'*weights(i)/2;
        rhs_elements2(:,i) = fespace.dettransf2* ...
            fespace.transffuns2{i}'*weights(i)/2;
    end
    for i = 1:n_elements
        indices = connectivity(i,1:end-1);
        x1 = vertices(indices(1),1:2)';
        
        if (indices(2) == indices(1) + 1)
            b(indices) = b(indices) + rhs_elements1* ...
                fun(fespace.transf1(gp,x1))';
            
        else
            b(indices) = b(indices) + rhs_elements2* ...
                fun(fespace.transf2(gp,x1))';
        end
    end
end