function [ Uti_max,Uti_min ] = Find_Uti_max_min( sigma,sigma_g,beta,B,r_d,N )
%FIND_UTI_MAX 此处显示有关此函数的摘要
%   此处显示详细说明
sigma_sum_max = max(sigma) * N;
sigma_sum_min = min(sigma) * N;
Uti = zeros((sigma_sum_max-sigma_sum_min)/(sigma(2)-sigma(1)),1);
i = 1;

for sum = sigma_sum_min: sigma(2)-sigma(1):sigma_sum_max
    Uti(i) = Utility_Reporter(sum,sigma_g,beta,B,r_d );
    i = i + 1; 
end

Uti_max = max(Uti);
Uti_min = min(Uti);
end

