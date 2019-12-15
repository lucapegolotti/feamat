function b = apply_dirichlet_bc_rhs(b,fespace,dirichlet_functions, varargin)
% Apply Dirichlet boundary conditions to rhs by evaluating the Dirichlet
% data in the Dirichlet nodes
%
% input=
%           b: right handside
%           fespace: finite element space
%           dirichlet_functions: boundary data
%           varargin: contains time, if unsteady problem is solved
%                                
% output=
%           b: right handside with boundary conditions

n_nodes = size(b,1);

nodes = fespace.nodes;
bc_flags = fespace.bc;

for i = 1:n_nodes
    if ((nodes(i,3)~=0) && bc_flags(nodes(i,3)))
        if nargin > 3
            vd = dirichlet_functions(nodes(i,1:2)', varargin{1});
        else
            vd = dirichlet_functions(nodes(i,1:2)');
        end
        b(i) = vd(nodes(i,3));
    elseif (nodes(i,4) ~= 0 && bc_flags(nodes(i,4)))
        if nargin > 3
            vd = dirichlet_functions(nodes(i,1:2)', varargin{1});
        else
            vd = dirichlet_functions(nodes(i,1:2)');
        end
        b(i) = vd(nodes(i,4));
    end
end

end
