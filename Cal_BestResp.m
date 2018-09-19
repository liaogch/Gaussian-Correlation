function [ BR, alpha ] = Cal_BestResp( N,j,sigma,cor_var,r,r_d )
%FIND_BESTRESP 此处显示有关此函数的摘要
%   此处显示详细说明
%M: number
%N: vector
%sigma: vector of all variance except sigma_j
%sum_cor_var: number
%r:social relation matrix
%

%sigma(j) = 0;
beta = 0;
alpha = 0;
for i=N
    beta = beta + r(j,i)*exp(-sum(sigma)- sum(cor_var(:,i)));
end

BR = max(0,-log(r_d/beta));

sigma_all = [sigma BR];
for i = N
    alpha = alpha + r(j,i)*exp(-sum(sigma_all)- sum(cor_var(:,i)));
end

