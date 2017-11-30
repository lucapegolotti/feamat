function [I] = compute_integral_over_mesh(u,fespace)
    v = assemble_vector_lagrange_multiplier(fespace);
    I = v'*u;
end

