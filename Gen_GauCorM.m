function [ L, W ] = Gen_GauCorM( T,Wmin,Wmax)
W = zeros(T);
D = zeros(T);

%mu = [1 2 3 2 1];
%sigma = [1 1.5 2 1.5 1];
%imax = [2 3 4];

for i = 1:T
    %W(i,:) = mu(i)+sigma(i)^2*randn(1,T);
    W(i,i+1:T) = Wmin+(Wmax-Wmin)*rand(1,T-i);
    W(i+1:T,i) = W(i,i+1:T)';
    W(i,i) = 0;
    D(i,i) = sum(W(i,:));
end

L = D-W;


end

