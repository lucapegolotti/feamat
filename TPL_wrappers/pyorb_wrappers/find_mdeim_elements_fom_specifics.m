function [elements] = find_mdeim_elements_fom_specifics( fem_specifics, indices_mat )
% Wrapper for extracting the list of mesh elements contributing to
% populating the given indices of the stiffness matrix, used in MDEIM
% framework
% input=
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           elements: array of elements

    [~, fespace] = set_fem_simulation( fem_specifics );

    elements = find_elements_given_matrix_indices( fespace, fespace, indices_mat );

end
