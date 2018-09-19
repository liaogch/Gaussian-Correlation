function [ res_sigma ] = BestResponseUpdate( beta,sigma_g,sigma )
%BESTRESPONSEUPDATE 此处显示有关此函数的摘要
%   此处显示详细说明
M  = length(sigma);
res_sigma = zeros(1,M);
for iter = 1:100
for i = 1:M
    res_sigma(i) = max(0,beta(i)-sigma_g-sum(sigma([1:i-1 i+1:M])));
end
sigma = res_sigma;
end


end

