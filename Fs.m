function [result] = Fs(H,theta,sigma,a,b,c)
% The function used to define the estimators
    h=0.1;
    fun = @(x) x.^(1-2*H)./(theta^2 + x.^2) ;
    fun1 = @(x) cos(h*x).*x.^(1-2*H)./(theta^2 + x.^2) ;
    fun2 = @(x) cos(10*h*x).*x.^(1-2*H)./(theta^2 + x.^2) ;
    f1 = sigma^2*gamma(2*H+1)*sin(pi*H)*integral(fun,0,100)/pi-a ;
    f2 = sigma^2*gamma(2*H+1)*sin(pi*H)*integral(fun1,0,100)/pi-b ;
    f3 = sigma^2*gamma(2*H+1)*sin(pi*H)*integral(fun2,0,100)/pi-c  ;
    result = [f1,f2,f3];
end

