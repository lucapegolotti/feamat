function [S,A,b,v,vt] = take_snapshots_poisson(fespace, fun, rangemu, nmus, rangeomega, nomegas, boundary_index)
% Take snapshots for the poisson problem where the parameter is the
% diffusion coefficient (constant over the domain). The additional
% parameter sets the integral on the non constrained boundary
% input= 
%           fespace: finite element space
%           rangemu: range of sampling of the diffusion parameter
%           nmus: number of samples taken from rangemu
%           rangeomega: range of sampling for omega (integral over non
%               constrained boundary)
%           nomegas: number of samples taken from rangeomega
%           boundary_index: index of the unconstrained boundary
% output=
%           S: matrix of snapshots

mus = linspace(rangemu(1),rangemu(2),nmus);
omegas = linspace(rangeomega(1),rangeomega(2),nomegas);

dirichlet_functions = @(x) [0;0;0;0];
neumann_functions = @(x) [0;0;0;0];
[v,vt] = assemble_vector_lagrange_multiplier_boundary(fespace,boundary_index);
S = zeros(size(fespace.nodes,1),nmus*nomegas);
count = 0;
for curmu = mus
    for curomega = omegas
        count = count + 1;
        disp(['Taking snapshot number ',num2str(count)]);
        [A,b] = assembler_poisson(fespace,fun,@(x) curmu,dirichlet_functions,neumann_functions);
        
        mat = [A v;vt 0];
        rhs = [b;curomega];
        
        sol = mat\rhs;
        S(:,count) = sol(1:end-1);
    end
end

A = A/curmu;

end

