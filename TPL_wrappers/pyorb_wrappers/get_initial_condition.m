function [u0] = get_initial_condition( param, fem_specifics, varargin)
% Get the initial condition of the problem, if the problem is unsteady
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the time marching scheme
%           varargin: test case number (optional)
% output=
%           u0: initial condition (if the problem is unsteady, NaN otherwise)

if (isfield(fem_specifics, 'final_time'))

    if nargin > 2
            caso = varargin{1};
        else
            caso = 1;
    end

    switch caso

        case 1
            %boundary conditions
            bc = [1;1;0;1];

            % initial condition
            u_init = @(x) 0*x(:,1) + 0*x(:,2);

        case 2
            %boundary conditions
            bc = [1;1;1;1];

            % initial condition
            u_init = @(x) (x(:,1)-x(:,1).^2).*(x(:,2)-x(:,2).^2);
    end
    [~, fespace] = set_fem_simulation( fem_specifics, bc );

    % evaluation of the initial condition
    u0 = u_init(fespace.nodes(:,1:2));

else
    
    u0 = nan; %initial condition is not defined if the problem is steady

end

           
end