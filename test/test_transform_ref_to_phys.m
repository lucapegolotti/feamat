clear all
close all
clc

x1 = [3;2];
x2 = [5;2.5];
x3 = [2.5;5];


A = [x2-x1 x3-x1];
b = x1;
detA = det(A);
f = @(x) A*x + b;

pt = [0.2;0.5];

subplot(1,2,1)
plot([0 1 0 0],[0 0 1 0]);
hold on
plot(pt(1),pt(2),'.','Markersize',20);
axis equal

subplot(1,2,2)
plot([x1(1) x2(1) x3(1) x1(1)], [x1(2) x2(2) x3(2) x1(2)])
trpt = f(pt);
hold on
plot(trpt(1),trpt(2),'.','Markersize',20);

axis equal

miny = @(x) (x<=x1(1)).*((x-x3(1))./(x1(1)-x3(1))*(x1(2)-x3(2)) + x3(2)) + ...
            (x>x1(1)).*((x-x1(1))./(x2(1)-x1(1))*(x2(2)-x1(2)) + x1(2));
        
maxy = @(x) ((x-x3(1))./(x2(1)-x3(1))*(x2(2)-x3(2)) + x3(2));
        
plot(2.5:0.01:5,miny(2.5:0.01:5));
plot(2.5:0.01:5,maxy(2.5:0.01:5));

fun = @(x,y) 1*x.^2 + y.*x;
% fun = @(x,y) sin(x).*y.^2;
funv = @(x) fun(x(1),x(2));
fun2 = @(x) funv(f(x));

integral2(fun,2.5,5,miny,maxy)

% integral2(fun,0,1,0,@(x) 1-x)

gp = [1/6 1/6; 2/3 1/6; 1/6 2/3];
weights = [1/3 1/3 1/3];

resint = 0;
for i = 1:3
    resint = resint + abs(det(A))*fun2(gp(i,:)') * weights(i) / 2;
end
resint


