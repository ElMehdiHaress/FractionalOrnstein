 Ltheta = zeros(10,1);
Ltheta1 = zeros(10,1);
Ltheta2 = zeros(10,1); 

Lsigma = zeros(10,1);
Lsigma1 = zeros(10,1);
Lsigma2 = zeros(10,1);

LH = zeros(10,1); 
LH1 = zeros(10,1);
LH2 = zeros(10,1);

theta0=4;
H0=0.6;
sigma0=3;

%{
 % >>> HEAVY BLOCK (commented out for CI) <<<
for s = 1:10
    N = Newton2(H0,theta0,sigma0,0.0001,20,s+6);
    n1 = length(N(1,:)) ;
    n2 = length(N(2,:)) ;
    n3 = length(N(3,:));
    theta = N(1,n1) ;
    sigma = N(2,n2);
    H = N(3,n3) ;
    Ltheta(s) = Ltheta(s) + theta ;
    LH(s) = LH(s) + H ;
    Lsigma(s) = Lsigma(s) + sigma ;
    
    N1 = Newton2(H0,theta0,sigma0,0.0001,20,s+6);
    n11 = length(N1(1,:)) ;
    n21 = length(N1(2,:)) ;
    n31 = length(N1(3,:)) ;
    theta1 = N1(1,n11) ;
    sigma1 = N1(2,n21);
    H1 = N1(3,n31) ;
    Ltheta1(s) = Ltheta1(s) + theta1 ;
    Lsigma1(s) = Lsigma1(s) + sigma ;
    LH1(s) = LH1(s) + H1 ;
    
    N2 = Newton2(H0,theta0,sigma0,0.0001,20,s+6);
    n12 = length(N2(1,:)) ;
    n22 = length(N2(2,:)) ;
    n32 = length(N2(3,:));
    theta2 = N2(1,n12) ;
    sigma2 = N2(2,n22);
    H2 = N2(3,n32) ;
    Ltheta2(s) = Ltheta2(s) + theta2 ;
    Lsigma2(s) = Lsigma2(s) + sigma2 ;
    LH2(s) = LH2(s) + H2 ;
        
end

plot(7:16,LH)
hold
plot(7:16,LH1)
plot(7:16,LH2)
%}

