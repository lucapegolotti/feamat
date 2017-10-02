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

elseif (strcmp('ROS3Pw low',label))
    gammad = 7.8867513459481287*1e-1;
    alpha = zeros(3);
    alpha(2,1) = 1.5773502691896257;
    alpha(3,1) = 0.5;
    
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
    coeffs.b = bhat;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('ROS2',label))
    gammad = 2.928932188*1e-1;
    alpha = zeros(2);
    alpha(2,1) = 1;

    b = zeros(2,1);
    b(1) = 0.5;
    b(2) = 0.5;

    bhat(1) = 0;
    bhat(2) = 1;

    gamma = zeros(2,2);
    gamma(2,1) = -0.58578643762;
    gamma(1,1) = gammad;
    gamma(2,2) = gammad;

    coeffs.alpha = alpha;
    coeffs.gamma = gamma;
    coeffs.b = b;
    coeffs.bhat = bhat;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('ROS3P',label))
    gammad =  7.886751346999999e-1;
    alpha = zeros(3);
    alpha(2,1) = 1;
    alpha(3,1) = 1;

    b = zeros(3,1);
    b(1) = 2/3;
    b(2) = 0;
    b(3) = 1/3;
    
    bhat = zeros(3,1);
    bhat(1) = 1/3;
    bhat(2) = 1/3;
    bhat(3) = 1/3;
    
    gamma = zeros(3,3);
    gamma(1,1) = gammad;
    gamma(2,2) = gammad;
    gamma(3,3) = gammad;
    gamma(2,1) = -1;
    gamma(3,1) = -7.886751346999999e-1;
    gamma(3,2) = -1.077350269000000;
    
    coeffs.alpha = alpha;
    coeffs.gamma = gamma;
    coeffs.b = b;
    coeffs.bhat = bhat;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('ROWDAIND2',label))
    gammad = 0.3;
    alpha = zeros(4);
    alpha(2,1) = 0.5;
    alpha(3,1) = 0.28;
    alpha(3,2) = 0.72;
    alpha(4,1) = 0.28;
    alpha(4,2) = 0.72;

    b = zeros(4,1);
    b(1) = 2/3;
    b(2) = 0;
    b(3) = 1/30;
    b(4) = 0.3;
    
    bhat = zeros(4,1);
    bhat(1) = 4.799002800355166e-1;
    bhat(2) = 5.176203811215082e-1;
    bhat(3) = 2.479338842975209e-3;
    bhat(4) = 0;
    
    gamma = zeros(4,4);
    gamma(1,1) = gammad;
    gamma(2,2) = gammad;
    gamma(3,3) = gammad;
    gamma(4,4) = gammad;

    gamma(2,1) = -1.121794871794876e-1;
    gamma(3,1) = 2.54;
    gamma(3,2) = -3.84;
    gamma(4,1) = 3.8666666666666666e-1;
    gamma(4,2) = -7.2e-1;
    gamma(4,3) = 1/30;
    
    coeffs.alpha = alpha;
    coeffs.gamma = gamma;
    coeffs.b = b;
    coeffs.bhat = bhat;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('ROS34PW2',label))
    gammad = 4.3586652150845900e-1;
    alpha = zeros(4);
    alpha(2,1) = 8.7173304301691801e-1;
    alpha(3,1) = 8.4457060015369423e-1;
    alpha(3,2) = -1.1299064236484185e-1;
    alpha(4,1) = 0;
    alpha(4,2) = 0;
    alpha(4,3) = 1;

    b = zeros(4,1);
    b(1) = 2.4212380706095346e-1;
    b(2) = -1.2232505839045147;
    b(3) = 1.5452602553351020;
    b(4) = 4.3586652150845900e-1;
    
    bhat = zeros(4,1);
    bhat(1) = 3.7810903145819369e-1;
    bhat(2) = -9.6042292212423178e-2;
    bhat(3) = 0.5;
    bhat(4) = 2.1793326075422950e-1;
    
    gamma = zeros(4,4);
    gamma(1,1) = gammad;
    gamma(2,2) = gammad;
    gamma(3,3) = gammad;
    gamma(4,4) = gammad;

    gamma(2,1) = -8.7173304301691801e-1;
    gamma(3,1) = -9.0338057013044082e-1;
    gamma(3,2) = 5.4180672388095326e-2;
    gamma(4,1) = 2.4212380706095346e-1;
    gamma(4,2) = -1.2232505839045147;
    gamma(4,3) = 5.4526025533510214e-1;
    
    coeffs.alpha = alpha;
    coeffs.gamma = gamma;
    coeffs.b = b;
    coeffs.bhat = bhat;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('ROS34PW3',label))
    gammad = 1.0685790213016289;
    alpha = zeros(4);
    alpha(2,1) = 2.5155456020628817;
    alpha(3,1) = 5.0777280103144085e-1;
    alpha(3,2) = 7.5000000000000000e-1;
    alpha(4,1) = 1.3959081404277204e-1;
    alpha(4,2) = -3.3111001065419338e-1;
    alpha(4,3) = 8.2040559712714178e-1;

    b = zeros(4,1);
    b(1) = 2.2047681286931747e-1;
    b(2) = 2.7828278331185935e-3;
    b(3) = 7.1844787635140066e-3;
    b(4) = 7.6955588053404989e-1;
    
    bhat = zeros(4,1);
    bhat(1) = 3.1300297285209688e-1;
    bhat(2) = -2.8946895245112692e-1;
    bhat(3) = 9.7646597959903003e-1;
    bhat(4) = 0;
    
    gamma = zeros(4,4);
    gamma(1,1) = gammad;
    gamma(2,2) = gammad;
    gamma(3,3) = gammad;
    gamma(4,4) = gammad;

    gamma(2,1) = -2.5155456020628817;
    gamma(3,1) = -8.7991339217106512e-1;
    gamma(3,2) = -9.6014187766190695e-1;
    gamma(4,1) = -4.1731389379448741e-1;
    gamma(4,2) = 4.1091047035857703e-1;
    gamma(4,3) = -1.3558873204765276;
    
    coeffs.alpha = alpha;
    coeffs.gamma = gamma;
    coeffs.b = b;
    coeffs.bhat = bhat;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('ROS3PRw',label))
    gammad = 7.8867513459481287e-1;
    alpha = zeros(3);
    alpha(2,1) = 2.3660254037844388;
    alpha(3,1) = 0.5;
    alpha(3,2) = 7.6794919243112270e-1;
    
    b = zeros(3,1);
    b(1) = 5.0544867840851759*1e-1;
    b(2) = -1.1571687603637559*1e-1;
    b(3) = 6.1026819762785800*1e-1;
    
    bhat = zeros(3,1);
    bhat(1) = -2.8973180237214197*1e-1;
    bhat(2) = 0.1;
    bhat(3) = 6.1026819762785800*1e-1;
    
    gamma = zeros(3,3);
    gamma(1,1) = gammad;
    gamma(2,2) = gammad;
    gamma(3,3) = gammad;
    gamma(2,1) = -2.3660254037844388;
    gamma(3,1) = -8.6791218280355165*1e-1;
    gamma(3,2) = -8.7306695894642317*1e-1;
    
    coeffs.alpha = alpha;
    coeffs.gamma = gamma;
    coeffs.b = b;
    coeffs.bhat = bhat;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('ROSI2P2',label))
    warning('ROSI2P2 seems not to converge witg the correct order')
    gammad = 4.3586652150845900e-1;
    alpha = zeros(4);
    alpha(2,1) = 0.5;
    alpha(3,1) = -5.1983699657507165e-1;
    alpha(3,2) = 1.5198369965750715;
    alpha(4,1) = -5.1983699657507165e-1;
    alpha(4,2) = 1.5198369965750715;

    b = zeros(4,1);
    b(1) = 2/3;
    b(2) = -5.4847955522165341e-32;
    b(3) = -1.0253318817512568e-1;
    b(4) = 4.3586652150845900e-1;
    
    bhat = zeros(4,1);
    bhat(1) = -9.5742384859111473e-1;
    bhat(2) =  2.9148476971822297;
    bhat(3) = 0.5;
    bhat(4) = -1.4574238485911146;
    
    gamma = zeros(4,4);
    gamma(1,1) = gammad;
    gamma(2,2) = gammad;
    gamma(3,3) = gammad;
    gamma(4,4) = gammad;

    gamma(2,1) = -0.5;
    gamma(3,1) = -4.0164172503011392e-1;
    gamma(3,2) = 1.1742718526976650;
    gamma(4,1) = 1.1865036632417383;
    gamma(4,2) = -1.5198369965750715;
    gamma(4,3) = -1.0253318817512568e-1;
    
    coeffs.alpha = alpha;
    coeffs.gamma = gamma;
    coeffs.b = b;
    coeffs.bhat = bhat;
    coeffs.nstages = size(alpha,1);
elseif (strcmp('ROSI2P1',label))
    gammad = 4.3586652150845900e-1;
    alpha = zeros(4);
    alpha(2,1) = 0.5;
    alpha(3,1) =  5.5729261836499822e-1;
    alpha(3,2) = 1.9270738163500176e-1;
    alpha(4,1) = -3.0084516445435860e-1;
    alpha(4,2) = 1.8995581939026787;
    alpha(4,3) = -5.9871302944832006e-1;

    b = zeros(4,1);
    b(1) =  5.2900072579103834e-2;
    b(2) =  1.3492662311920438;
    b(3) =  -9.1013275270050265e-1;
    b(4) =  5.0796644892935516e-1;
    
    bhat = zeros(4,1);
    bhat(1) =  1.4974465479289098e-1;
    bhat(2) =  7.0051069041421810e-1;
    bhat(3) = 0;
    bhat(4) =  1.4974465479289098e-1;
    
    gamma = zeros(4,4);
    gamma(1,1) = gammad;
    gamma(2,2) = gammad;
    gamma(3,3) = gammad;
    gamma(4,4) = gammad;

    gamma(2,1) = -0.5;
    gamma(3,1) = -6.4492162993321323e-1;
    gamma(3,2) = 6.3491801247597734e-2;
    gamma(4,1) = 9.3606009252719842e-3;
    gamma(4,2) = -2.5462058718013519e-1;
    gamma(4,3) = -3.2645441930944352e-1;
    
    coeffs.alpha = alpha;
    coeffs.gamma = gamma;
    coeffs.b = b;
    coeffs.bhat = bhat;
    coeffs.nstages = size(alpha,1);
else
    error('Rosenbrock method not implemented');
end