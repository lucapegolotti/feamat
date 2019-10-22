function b_adj = adjust_rhs(b,f_t,fespace)
% Function to adjust the rhs vector of the parabolic problem, in order to
% take into account the presence of the time dependent part of the forcing
% term
%
% Input=
%            b: original rhs vector, with bc already imposed
%            f_t: function handle, containing the time dependent expresion of
%            the forcing term
%            fespace: finite element space structure
%
% Output=
%              b_adj: function handle (in time) containing the adjusted 
%              evaluation of the rhs vector 

n_nodes = size(b,1);
nodes = fespace.nodes;

for i = 1:n_nodes
    if (nodes(i,3)~=0)
        b_adj= @(t) b / f_t(t);
    end
end

end