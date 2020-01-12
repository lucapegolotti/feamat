function [] = plot_fe_solution(rb_solution, fe_solution, fem_specifics, folder, name, varargin)
% Plot the RB solution and the FE solution of a problem
% input=
%           rb_solution vector containing the RB solution to the problem at
%           hand
%           fe_solution: vector containing the FE solution to the problem at
%           hand. If empty vector or double, just the RB solution is plotted
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the time marching scheme
%           folder: folder where to save the plot
%           name: name of the plot
%           varargin: name of the test case, if specified. Else it is set
%           to 1

% dbstop at 32

if nargin > 3
        caso = varargin{1};
else
        caso = 1;
end

switch caso
    case 1
        bc_flags = [1;1;0;1];
    case 2
        bc_flags = [1;1;1;1];
end

[~, fespace] = set_fem_simulation(fem_specifics, bc_flags);

fig = figure();
if ~isempty(fe_solution)
    for i = 1:2
        subplot(1,2,i)
        if i == 1
            plot_fe_function(rb_solution, fespace);
            title('RB Solution Plot')
        elseif i==2
            plot_fe_function(fe_solution, fespace);
            title('FE Solution Plot')
        end
    end
else
   plot_fe_function(rb_solution, fespace);
   title('Error Plot') 
end

folder= convertCharsToStrings(folder);
name = convertCharsToStrings(name);
             
saveas(fig, folder+name, 'epsc')
%savefig(folder+name+'.fig')

end