function [ s ] = Gen_SocRelMat( T,M,N,Smin,Smax )

s= zeros(T);

ss = Smin+(Smax-Smin)*rand(length(M),length(N));

s(M,N)=ss;

end

