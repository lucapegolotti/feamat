function draw_mesh_element(index,vertices,connectivity)

x = [vertices(connectivity(index,1),1) vertices(connectivity(index,2),1) vertices(connectivity(index,3),1)]; 
y = [vertices(connectivity(index,1),2) vertices(connectivity(index,2),2) vertices(connectivity(index,3),2)];

plot([x x(1)],[y y(1)],'k');