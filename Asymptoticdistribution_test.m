Ltheta = zeros(100,1);
 

Lsigma = zeros(100,1);


LH = zeros(100,1); 


theta0=4;
H0=0.6;
sigma0=3;

%{
  % >>> HEAVY BLOCK (commented out for CI) <<<
for s = 1:100
    N = Newton2(H0,theta0,sigma0,0.0001,20,12);
    n1 = length(N(1,:)) ;
    n2 = length(N(2,:)) ;
    n3 = length(N(3,:));
    theta = N(1,n1) ;
    sigma = N(2,n2);
    H = N(3,n3) ;
    Ltheta(s) = Ltheta(s) + theta ;
    LH(s) = LH(s) + H ;
    Lsigma(s) = Lsigma(s) + sigma ;
    
end

[f,xi] = ksdensity(Ltheta); 
figure
hold
plot(xi,f);

[f,xi] = ksdensity(Lsigma); 
figure
plot(xi,f);

[f,xi] = ksdensity(LH); 
figure
plot(xi,f)

M1 = length(Ltheta);
M2 = length(Lsigma);
M3 = length(LH);

X = sum(Ltheta)/M1 ;
Y = sum(Lsigma)/M2 ;
S = sum(LH)/M3 ;

Vx = sum((Ltheta-X).^2)/(M1-1) ;
Vy = sum((Lsigma-Y).^2)/(M2-1) ;
Vs = sum((LH-S).^2)/(M3-1) ;

disp([[X,Y,S];[sqrt(Vx),sqrt(Vy),sqrt(Vs)]])
%}


