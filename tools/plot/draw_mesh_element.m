function draw_mesh_element(index,vertices,elements,varargin)
% Draw a single mesh element
% input=
%           index: index of the element
%           vertices: vertices of the element
%           elements: elements of the mesh
%           (optional)
%           'color' to color the boundary elements with different colors

group = elements(index,4);

x = [vertices(elements(index,1),1) vertices(elements(index,2),1) vertices(elements(index,3),1)]; 
y = [vertices(elements(index,1),2) vertices(elements(index,2),2) vertices(elements(index,3),2)];

color = 0;
if (nargin == 4)
    if (strcmp(varargin{1},'color'))
        color = 1;
    else
        error('Unknown parameter to draw_mesh_element function!');
    end
end

if (~color)
    plot([x x(1)],[y y(1)],'color',[0.7 0.7 0.7],'Linewidth',0.5);
else
    if (group == 1)
        plot([x x(1)],[y y(1)],'r','Linewidth',1);
    elseif (group == 2)
        plot([x x(1)],[y y(1)],'b','Linewidth',1);
    elseif (group == 3)
        plot([x x(1)],[y y(1)],'g','Linewidth',1);
    elseif (group == 4)
        plot([x x(1)],[y y(1)],'y','Linewidth',1);
    elseif (group == 0)
        plot([x x(1)],[y y(1)],'k','Linewidth',0.5)
    end
end
