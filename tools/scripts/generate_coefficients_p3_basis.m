clear all
clc

% compute coefficients of basis functions for P3 polynomials


row = @(x) [x(1)^3 x(2)^3 x(1)^2*x(2) x(1)*x(2)^2 x(1)^2 x(2)^2 x(1)*x(2) x(1) x(2) 1];

a = 1/3;
b = 2/3;
p1 = [0 0];
p2 = [1 0];
p3 = [0 1];
p4 = [a 0];
p5 = [b 0];
p6 = [b a];
p7 = [a b];
p8 = [0 b];
p9 = [0 a];
p10 = (p1 + p2 + p3)/3;

A = [row(p1); row(p2); row(p3); row(p4); row(p5); row(p6); ...
     row(p7); row(p8); row(p9); row(p10)];


m = zeros(10);
string_matrix = 'c = [';
for i = 1:10
    b = zeros(10,1);
    b(i) = 1;
    m(i,:) = A\b;
    for j = 1:10
        if (abs(m(i,j)) < 1e-10)
            %m(i,j) = 0;
        end
        string_matrix = [string_matrix,num2str(m(i,j),'%10.5e\n')];
        if (j~=10)
            string_matrix = [string_matrix,','];
        end
    end
    if (i ~= 10)
        string_matrix = [string_matrix,';'];
    end
end
string_matrix = [string_matrix,'];'];
disp(string_matrix)