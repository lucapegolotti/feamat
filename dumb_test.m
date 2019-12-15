clear vars
close all
clc
%%

params = [1.00, 2.00, 3.00];
spatial_dim = 50;

fem_specifics.number_of_elements = spatial_dim;
fem_specifics.polynomial_degree = 'P1';
fem_specifics.model = 'thermal_block';
fem_specifics.use_nonhomogeneous_dirichlet = 'N';                                                            
fem_specifics.mesh_name = ['cube' num2str(spatial_dim) 'x' num2str(spatial_dim)];

caso = 1;

% New fields in time dependent case!!!!
switch caso
    case 1
        fem_specifics.final_time = 0.1;
    case 2
        fem_specifics.final_time = 1.0;
end
fem_specifics.number_of_time_instances = 50;
fem_specifics.theta = 1.0;

timestep_number = 4;
test_nb = 1;

sol = get_reference_sol(params, fem_specifics, timestep_number, caso);
