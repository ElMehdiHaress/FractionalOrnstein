function [ result ] = Newton2( H0,theta0,sigma0, tol, nmax,s )
%
% NEWTON Newton's Method
%   Newton's method for finding successively better approximations to the 
%   zeroes of a real-valued function.
%
% Input:
%   H0, theta0, sigma0 - true parameters
%   tol - tolerance
%   nmax - maximum number of iterations
%    s  -  log_2 of the sample size
% Output:
%   theta, H, sigma - estimates
%   ex - error estimate
%
%


    
    
    h = 0.1 ;
    a = eta_quantities(s,h) ;
    b = eta_quantities(s,10*h) ;
    a1=a(1);
    a2=a(2);
    a3=b(2);
    
    F = Fs(H0,theta0,sigma0,a1,a2,a3);
    f1 = F(1);
    f2 = F(2);
    f3 = F(3) ;
   
    D = Determinant(theta0,H0,sigma0) ;
    Delta = linsolve(D,-[f1;f2;f3]);  
    
    theta(1) = real(theta0 + Delta(1)) ;
    H(1) = abs(H0+Delta(3)) ;
    sigma(1) = abs(sigma0+Delta(2)) ;
    ex(1) = abs((theta(1)-theta0)^2 + (sigma(1)-sigma0)^2 + (H(1)-H0)^2);

    k = 2;
    while  (ex(k-1)>= tol) && (k <= nmax) %(ex(k-1) >= tol) &&
%         while Dfunctionf(H(k-1)) <= 0
%         epsilon = rand()*0.1 - 0.1 ;
%         H(k-1) = H(k-1) + epsilon ;
%         end
    F = Fs(H(k-1),theta(k-1),sigma(k-1),a1,a2,a3);
    f1 = F(1);
    f2 = F(2);
    f3 = F(3) ;
   
    D = Determinant(theta(k-1),H(k-1),sigma(k-1)) ;
    Delta = linsolve(D,-[f1;f2;f3]);     
    
    theta(k) = real(theta(k-1) + Delta(1)) ;
    H(k) = abs(H(k-1)+Delta(3)) ;
    sigma(k) = abs(sigma(k-1)+Delta(2)) ;
    ex(k) = abs((theta(k)-theta(k-1))^2 + (sigma(k)-sigma(k-1))^2 + (H(k)-H(k-1))^2);
        k = k+1;
    end

result = [theta;sigma;H;ex] ;

end
