function [fespace, gp, weights] = add_tensors_structured_meshes(fespace, n_gauss)
    % Add members to fespace which are necessary to optimize assembly on
    % structured meshes. Such members are evaluations of terms of the type
    % grad(phi_j)_p*grad(phi_i)_q, i.e. the nondimensional (p,q) components
    % of the diffusivity tensor (p,q=1,2).
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
    fespace.diffusivity_elements1 = cell(2,2);
    fespace.diffusivity_elements_sum1 = cell(2,2);
    invmatT = (fespace.mattransf1 \ eye(2))';

    % assemble (p,q) components of the diffusivity tensor separately
    % only the upper triangular part is saved
    for i = 1:n_gauss
        transformedGrad = invmatT * fespace.grads(gp(:,i));
        for p = 1:2
            for q = 1:2
                fespace.diffusivity_elements1{p,q}{end+1} = fespace.dettransf1 * ...
                    transformedGrad(p,:)' * transformedGrad(q,:) * weights(i)/2;
            end
        end
    end

    for p = 1:2
        for q = 1:2
            fespace.diffusivity_elements_sum1{p,q} = fespace.diffusivity_elements1{p,q}{1};
            for i = 2:n_gauss
                fespace.diffusivity_elements_sum1{p,q} = fespace.diffusivity_elements_sum1{p,q} + ...
                    fespace.diffusivity_elements1{p,q}{i};
            end
        end
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
    fespace.diffusivity_elements2 = cell(2,2);
    fespace.diffusivity_elements_sum2 = cell(2,2);
    invmatT = (fespace.mattransf2 \ eye(2))';

    for i = 1:n_gauss
        transformedGrad = invmatT * fespace.grads(gp(:,i));
        for p = 1:2
            for q = 1:2
                fespace.diffusivity_elements2{p,q}{end+1} = fespace.dettransf2 * ...
                    transformedGrad(p,:)' * transformedGrad(q,:) * weights(i)/2;
            end
        end
    end

    for p = 1:2
        for q = 1:2
            fespace.diffusivity_elements_sum2{p,q} = fespace.diffusivity_elements2{p,q}{1};
            for i = 2:n_gauss
                fespace.diffusivity_elements_sum2{p,q} = fespace.diffusivity_elements_sum2{p,q} + ...
                    fespace.diffusivity_elements2{p,q}{i};
            end
        end
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