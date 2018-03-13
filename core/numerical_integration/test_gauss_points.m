% common variables

tol = 1e-12;

%% Test 1: 1D 2 points

[gp,weights,order] = gauss_points1D(2);

assert(order == 3);

coefs = rand(4,1);

f = @(x) coefs(1)*x.^0+coefs(2)*x.^1+coefs(3)*x.^2+coefs(4)*x.^3;
appr_int = weights*f(gp');

exact_int = integral(f,-1,1);
assert(abs(appr_int-exact_int) < tol);

%% Test 2: 1D 3 points

[gp,weights,order] = gauss_points1D(3);

assert(order == 5);

coefs = rand(6,1);

f = @(x) coefs(1)*x.^0+coefs(2)*x.^1+coefs(3)*x.^2+ ...
         coefs(4)*x.^3+coefs(5)*x.^4+coefs(6)*x.^5;
appr_int = weights*f(gp');

exact_int = integral(f,-1,1);
assert(abs(appr_int-exact_int) < tol);

%% Test 3: 1D 4 points

[gp,weights,order] = gauss_points1D(4);

assert(order == 7);

coefs = rand(8,1);

f = @(x) coefs(1)*x.^0+coefs(2)*x.^1+coefs(3)*x.^2+ ...
         coefs(4)*x.^3+coefs(5)*x.^4+coefs(6)*x.^5+ ...
         coefs(7)*x.^6+coefs(8)*x.^7;
appr_int = weights*f(gp');

exact_int = integral(f,-1,1);
assert(abs(appr_int-exact_int) < tol);

%% Test 4: 2D 3 points

[gp,weights,order] = gauss_points2D(3);

assert(order == 2);

f = @(x) x(:,1).^2 + x(:,1).*x(:,2) - 1;

appr_int = weights*f(gp')/2;

assert(abs(appr_int + 0.375) < tol)

%% Test 4: 2D 4 points

[gp,weights,order] = gauss_points2D(4);

assert(order == 3);

f = @(x) (x(:,1) + x(:,2)).^3;

appr_int = weights*f(gp')/2;

assert(abs(appr_int - 0.2) < tol)