function [I] = compute_integral_over_boundary(u,fespace,boundary_index)
    v = assemble_vector_lagrange_multiplier_boundary(fespace,boundary_index);
    I = v'*u;
end

