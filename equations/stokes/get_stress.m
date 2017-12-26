function [x,stress] = get_stress(sol,coor,label,epsilon)

[x,u1] = get_values_over_line(sol.fespace_u,sol.u1,100,coor,label);
[~,u2] = get_values_over_line(sol.fespace_u,sol.u2,100,coor,label);

[~,p] = get_values_over_line(sol.fespace_p,sol.p,100,coor,label);

[~,u1eps] = get_values_over_line(sol.fespace_u,sol.u1,100,coor+epsilon,label);
[~,u2eps] = get_values_over_line(sol.fespace_u,sol.u2,100,coor+epsilon,label);

du1 = (u1eps-u1)/epsilon;
du2 = (u2eps-u2)/epsilon;

if (strcmp(label,'Xpar'))
    stress = [du1 du2-p];
else
    stress = [du1-p du2];
end
