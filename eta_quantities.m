function [result] = eta_quantities(s,h) 
% This function computes the quantities needed to define the estimator


tic
n = 2^s;
N = 2^9 ;

% Insert the true parameters here:
H = 0.4 ;
Theta = 6 ;
sigma = 2 ;

B = fbm1d(H,N*n,n*h) ;


X = zeros(n+1,1) ;


for k = 2:N*n+1
    X(k) = exp(-Theta*h/N)*X(k-1) + sigma*exp(-Theta*h/N)*(B(k)-B(k-1)) ;
end



H2 = 0 ;
for k=1:n
       H2 = H2 + (X(k*N)^2)/n ;
       
end


H3 = 0 ;


for k = 1:(n-1)
       H3 = H3 + (X((k + 1)*N)*X(k*N))/n ;
end


    



result = [H2,H3];

toc

end












