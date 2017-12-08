clear all
clc

% compute coefficients of basis functions for P2 polynomials


row = @(x) [x(1)^2 x(2)^2 x(1)*x(2) x(1) x(2) 1];

p1 = [0 0];
p2 = [0.5 0];
p3 = [1 0];
p4 = [0.5 0.5];
p5 = [0 1];
p6 = [0 0.5];

A = [row(p1); row(p2); row(p3); row(p4); row(p5); row(p6)];


% coefficients for first basis function
b = zeros(6,1);
b(1) = 1;
A\b

% coefficients for second basis function
b = zeros(6,1);
b(2) = 1;
A\b

% coefficients for third basis function
b = zeros(6,1);
b(3) = 1;
A\b

% coefficients for fourth basis function
b = zeros(6,1);
b(4) = 1;
A\b

% coefficients for fifth basis function
b = zeros(6,1);
b(5) = 1;
A\b

% coefficients for sixth basis function
b = zeros(6,1);
b(6) = 1;
A\b
