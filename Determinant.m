function [result] = Determinant(theta,H,sigma)
% Jacobian matrix of Fs
    h = 0.1 ;
    fun1 = @(x) x.^(1-2*H)./ (theta^2+x.^2); 
    fun11 = @(x) -2.*log(x).*x.^(1-2*H)./ (theta^2+x.^2);
    fun12 = @(x) -2*theta.*x.^(1-2*H)./ (theta^2+x.^2).^2;
    fun13 = @(x) 2.*log(x).*x.^(2*H).*exp(-x) ;
    
    
    L1 = [(1/pi)*sigma^2 *gamma(2*H+1)*sin(pi*H)*integral(fun12,0,100), (2/pi)*sigma* gamma(2*H+1)*sin(pi*H)*integral(fun1,0,100), (1/pi)*sigma^2*(sin(pi*H)*integral(fun13,0,100)*integral(fun1,0,100) + pi*gamma(2*H+1)*cos(pi*H)*integral(fun1,0,100) + gamma(2*H+1)*sin(pi*H)*integral(fun11,0,100) ) ];
    
    fun2 = @(x) cos(h.*x).*x.^(1-2*H) ./ (theta^2+x.^2) ;
    fun21 = @(x) -2.*cos(h.*x).*log(x).*x.^(1-2*H)./(theta^2+x.^2) ;
    fun22 = @(x)  -2*theta.*cos(h.*x).*x.^(1-2*H) ./ (theta^2+x.^2).^2 ;
    fun23 = @(x) 2.*log(x).*x.^(2*H).*exp(-x) ;
    
    L2 = [(1/pi)*sigma^2 * gamma(2*H+1)*sin(pi*H)*integral(fun22,0,100),(2/pi)*sigma* gamma(2*H+1)*sin(pi*H)*integral(fun2,0,100), (1/pi)*sigma^2*(sin(pi*H)*integral(fun23,0,100)*integral(fun2,0,100) + pi*gamma(2*H+1)*cos(pi*H)*integral(fun2,0,100) + gamma(2*H+1)*sin(pi*H)*integral(fun21,0,100) ) ] ;
        
    fun3 = @(x) cos(10*h.*x).*x.^(1-2*H) ./ (theta^2+x.^2) ;
    fun31 = @(x) -2.*cos(10*h.*x).*log(x).*x.^(1-2*H)./(theta^2+x.^2) ;
    fun32 = @(x)  -2*theta.*cos(10*h.*x).*x.^(1-2*H) ./ (theta^2+x.^2).^2 ;
    fun33 = @(x) 2.*log(x).*x.^(2*H).*exp(-x) ;
    L3 = [(1/pi)*sigma^2 * gamma(2*H+1)*sin(pi*H)*integral(fun32,0,100),(2/pi)*sigma* gamma(2*H+1)*sin(pi*H)*integral(fun2,0,100), (1/pi)*sigma^2*(integral(fun33,0,100)*integral(fun3,0,100) + pi*gamma(2*H+1)*cos(pi*H)*integral(fun3,0,100) + gamma(2*H+1)*sin(pi*H)*integral(fun31,0,100) ) ] ;
    
    result = [L1;L2;L3] ;

end
    
    
    
