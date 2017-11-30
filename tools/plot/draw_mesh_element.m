function draw_mesh_element(index,vertices,elements,varargin)

group = elements(index,4);

x = [vertices(elements(index,1),1) vertices(elements(index,2),1) vertices(elements(index,3),1)]; 
y = [vertices(elements(index,1),2) vertices(elements(index,2),2) vertices(elements(index,3),2)];

if (nargin == 3)
    if (group == 1)
        plot([x x(1)],[y y(1)],'r','Linewidth',1);
    elseif (group == 2)
        plot([x x(1)],[y y(1)],'b','Linewidth',1);
    elseif (group == 3)
        plot([x x(1)],[y y(1)],'g','Linewidth',1);
    elseif (group == 4)
        plot([x x(1)],[y y(1)],'y','Linewidth',1);
    elseif (group == 0)
        plot([x x(1)],[y y(1)],'k');
    end
elseif (strcmp(varargin{1},'nocolor'))
    plot([x x(1)],[y y(1)],'k','Linewidth',1);
else
    error('Option is not supported');
end
