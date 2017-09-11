function coeffs = rosenbrock_coeffs(label)

if (strcmp('ROS3Pw',label))
    gammad = 7.8867513459481287*1e-1;
    alpha = zeros(3);
    alpha(2,1) = 1.5773502691896257;
    alpha(3,1) = 0.5;
    
    b = zeros(3,1);
    b(1) = 1.0566243270259355*1e-1;
    b(2) = 4.9038105676657971*1e-2;
    b(3) = 8.4529946162074843*1e-1;
    
    bhat = zeros(3,1);
    bhat(1) = -1.7863279495408180*1e-1;
    bhat(2) = 1/3;
    bhat(3) = 8.4529946162074843*1e-1;
    
    gamma = zeros(3,3);
    gamma(1,1) = gammad;
    gamma(2,2) = gammad;
    gamma(3,3) = gammad;
    gamma(2,1) = -1.5773502691896257;
    gamma(3,1) = -6.7075317547305480*1e-1;
    gamma(3,2) = -1.7075317547305482*1e-1;
    
    coeffs.alpha = alpha;
    coeffs.gamma = gamma;
    coeffs.b = b;
    coeffs.bhat = bhat;
    coeffs.nstages = size(alpha,1);
else
    error('Rosenbrock method not implemented');
end