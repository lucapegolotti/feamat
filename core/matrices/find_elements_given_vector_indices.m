function found_elements = find_elements_given_vector_indices( fespace, indices )

elements = fespace.connectivity;

found_elements = [];
for i = 1:length(indices)

    aux = elements(:,1) == indices(i);
    for j = 2:fespace.n_functions_per_element
        aux = aux | elements(:,j) == indices(i);
    end
    
    el = find(aux);
    found_elements = [found_elements; el];
end

% delete repetitions
found_elements = unique(found_elements);
