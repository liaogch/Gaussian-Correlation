function [ Uti ] = Utility_Reporter( sigma_sum,sigma_g,beta,B,r_d )
%UTILITY_REPORTER 此处显示有关此函数的摘要
%   此处显示详细说明

Uti = B - r_d * (sigma_sum+sigma_g) - exp(-sigma_sum-sigma_g)*beta;

end

