function elements = find_elements_given_matrix_indices(fespace_1,fespace_2,indices)

elements_1 = fespace_1.connectivity;
elements_2 = fespace_2.connectivity;

elements = [];
for i = 1:size(indices,1)

    aux = elements_1(:,1) == indices(i,1);
    for j = 2:fespace_1.n_functions_per_element
        aux = aux | elements_1(:,j) == indices(i,1);
    end
    
    el1 = find(aux);
    
    aux = elements_2(:,1) == indices(i,2);
    for j = 2:fespace_2.n_functions_per_element
        aux = aux | elements_2(:,j) == indices(i,2);
    end
    
    el2 = find(aux);
    
    new_elements = intersect(el1,el2);
    elements = [elements;new_elements];
end

% delete repetitions
elements = unique(elements);

length(elements)
