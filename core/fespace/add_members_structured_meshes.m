function [fespace,gp,weights] = add_members_structured_meshes(fespace,n_gauss)
% Add members to fespace which are necessary to optimize assembly on
% structured meshes
% input=
%           fespace: finite element space to modify
%           n_gauss: number of gauss points
%
% output=
%           fespace: modified finite element space


[gp,weights,~] = gauss_points2D(n_gauss);

% add structures for configuration 1, i.e. triangles with this orientation
%     /|
%    / |
%   /  |
%   ----

fespace.mattransf1 = [fespace.mesh.h1 fespace.mesh.h1; 0 fespace.mesh.h2];
fespace.dettransf1 = abs(det(fespace.mattransf1));
fespace.transf1 = @(x,x1) fespace.mattransf1 * x + x1;
fespace.transfgrads1 = {};
fespace.stiffness_elements1 = {};
invmat = inv(fespace.mattransf1);

for i = 1:n_gauss
    fespace.transfgrads1{end+1} = invmat' * fespace.grads(gp(:,i));
    fespace.stiffness_elements1{end+1} = fespace.dettransf1* ...
        (fespace.transfgrads1{i}'* ...
        fespace.transfgrads1{i})* ...
        weights(i)/2;
end

fespace.stiffness_elements_sum1 = fespace.stiffness_elements1{1};
for i = 2:n_gauss
    fespace.stiffness_elements_sum1 = fespace.stiffness_elements_sum1 + ...
        fespace.stiffness_elements1{i};
end

fespace.transffuns1 = {};
fespace.mass_elements1 = {};
for i = 1:n_gauss
    fespace.transffuns1{end+1} = fespace.functions(gp(:,i))';
    fespace.mass_elements1{end+1} = fespace.dettransf1* ...
        (fespace.transffuns1{i}'* ...
        fespace.transffuns1{i})* ...
        weights(i)/2;
end

fespace.mass_elements_sum1 = fespace.mass_elements1{1};
for i = 2:n_gauss
    fespace.mass_elements_sum1 = fespace.mass_elements_sum1 + ...
        fespace.mass_elements1{i};
end


% add structures for configuration 2, i.e. triangles with this orientation
%   ____
%   |  /
%   | /
%   |/

fespace.mattransf2 = [0 fespace.mesh.h1; fespace.mesh.h2 fespace.mesh.h2];
fespace.dettransf2 = abs(det(fespace.mattransf2));
fespace.transf2 = @(x,x1) fespace.mattransf2 * x + x1;
fespace.transfgrads2 = {};
fespace.stiffness_elements2 = {};
invmat = inv(fespace.mattransf2);

for i = 1:n_gauss
    fespace.transfgrads2{end+1} = invmat' * fespace.grads(gp(:,i));
    fespace.stiffness_elements2{end+1} = fespace.dettransf2* ...
        (fespace.transfgrads2{i}'* ...
        fespace.transfgrads2{i})* ...
        weights(i)/2;
end

fespace.stiffness_elements_sum2 = fespace.stiffness_elements2{1};
for i = 2:n_gauss
    fespace.stiffness_elements_sum2 = fespace.stiffness_elements_sum2 + ...
        fespace.stiffness_elements2{i};
end

fespace.transffuns2 = {};
fespace.mass_elements2 = {};
for i = 1:n_gauss
    fespace.transffuns2{end+1} = fespace.functions(gp(:,i))';
    fespace.mass_elements2{end+1} = fespace.dettransf2* ...
        (fespace.transffuns2{i}'* ...
        fespace.transffuns2{i})* ...
        weights(i)/2;
end

fespace.mass_elements_sum2 = fespace.mass_elements2{1};
for i = 2:n_gauss
    fespace.mass_elements_sum2 = fespace.mass_elements_sum2 + ...
        fespace.mass_elements2{i};
end

end