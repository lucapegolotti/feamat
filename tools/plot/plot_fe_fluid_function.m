function plot_fe_fluid_function(sol,what,varargin)
% Plot finite element fluid function
%
% input=
%           sol: fluid datastructured (created with solve_fluid_system)
%           what: 'U' or 'P', it specifies which field to plot (either
%           velocity or pressure)
%           (optional)
%           a 2 component vector that specifies the minimum and maximum
%           value to plot
%

if (nargin >= 3)
    arg = varargin{1};
    mlim = arg(1);
    Mlim = arg(2);
end

if (strcmp(sol.fespace_u.mesh.type,'structured'))
    if (what == 'U')

        fespace_u = sol.fespace_u;

        n_vertices = size(fespace_u.mesh.vertices,1);
        n1 = size(fespace_u.mesh.X,1);
        n2 = size(fespace_u.mesh.X,2);
        U1 = reshape(sol.u1(1:n_vertices),n1,n2);
        U2 = reshape(sol.u2(1:n_vertices),n1,n2);

        N = sqrt(U1.^2+U2.^2);

        [~,c] = contourf(fespace_u.mesh.X,fespace_u.mesh.Y,N);
        c.LineStyle = 'none';
        shading interp
        h = colorbar;
        if (nargin >= 3)
            if (mlim == 0 && Mlim == 0)
                mlim = 0; 
                Mlim = 1;
            end
            set(h,'YLim',[mlim*0.9 Mlim*1.1])
            caxis([ mlim*0.9 Mlim*1.1 ])
            c.LevelList = linspace(mlim,Mlim,20);
        else
            m = min(min(N));
            M = max(max(N));
            c.LevelList = linspace(m,M,20);
        end
        hold on

        q = quiver(fespace_u.mesh.X,fespace_u.mesh.Y,U1./N,U2./N);
        q.Color = 'black';
        q.LineWidth = 1;
        q.AutoScaleFactor = 0.5;

    elseif (what == 'P')
        fespace_p = sol.fespace_p;

        p = sol.p;

        n1 = size(fespace_p.mesh.X,1);
        n2 = size(fespace_p.mesh.X,2);

        [~,c] = contourf(fespace_p.mesh.X,fespace_p.mesh.Y,reshape(p,n1,n2));
        c.LineStyle = 'none';

        h = colorbar;  

        if (nargin >= 5)
            if (mlim == 0 && Mlim == 0)
                mlim = 0; 
                Mlim = 1;
            end
            set(h,'YLim',[mlim*0.9 Mlim*1.1])
            caxis([ mlim*0.9 Mlim*1.1 ])
        else
            m = min(p);
            M = max(p);
            c.LevelList = linspace(m,M,20);
        end
    end
elseif (strcmp(sol.fespace_u.mesh.type,'unstructured'))
    if (what == 'U')
        fespace_u = sol.fespace_u;

        n_vertices = size(fespace_u.mesh.vertices,1);
        N = sqrt(sol.u1(1:n_vertices).^2+sol.u2(1:n_vertices).^2);

        [~,h] = tricontf(fespace_u.mesh.vertices(:,1),fespace_u.mesh.vertices(:,2), ...
               fespace_u.mesh.elements(:,1:3),N,20);
        hold on
        %plot_fe_function(N,fespace_u);
        % set(h,'edgecolor','none');
%         h = triquiver(fespace_u.mesh.elements(:,1:3),fespace_u.mesh.vertices(:,1), ...
%                   fespace_u.mesh.vertices(:,2), sol.u1(1:n_vertices),sol.u2(1:n_vertices),10);
%         shading interp
%         h = colorbar;
%         if (nargin >= 3)
%             if (mlim == 0 && Mlim == 0)
%                 mlim = 0; 
%                 Mlim = 1;
%             end
%             set(h,'YLim',[mlim*0.9 Mlim*1.1])
%             caxis([ mlim*0.9 Mlim*1.1 ])
%             c.LevelList = linspace(mlim,Mlim,20);
%         else
%             m = min(min(N));
%             M = max(max(N));
%             c.LevelList = linspace(m,M,20);
%         end
%         hold on

%         q = quiver(fespace_u.mesh.X,fespace_u.mesh.Y,U1./N,U2./N);
%         q.Color = 'black';
%         q.LineWidth = 1;
%         q.AutoScaleFactor = 0.5;

    elseif (what == 'P')
        fespace_p = sol.fespace_p;
        
        [~,h] = tricontf(fespace_p.mesh.vertices(:,1),fespace_p.mesh.vertices(:,2), ...
               fespace_p.mesh.elements(:,1:3),sol.p,20);
        hold on
        set(h,'edgecolor','none');
        
    end
else
    error('Type of mesh is not supported!');
end