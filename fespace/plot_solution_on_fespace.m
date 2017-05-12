function plot_solution_on_fespace(fespace,vec)

n_vertices = size(fespace.mesh.vertices,1);
n1 = size(fespace.mesh.X,1);
n2 = size(fespace.mesh.X,2);
surf(fespace.mesh.X,fespace.mesh.Y,reshape(vec(1:n_vertices),n1,n2),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
