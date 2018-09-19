function [ Ave_Uti  ] = Cal_AveUti_Rep( M,T, sum_sigma_d,sigma_g,cor_var,S,R_d,r_d )
%CAL_AVEUTI_REP 此处显示有关此函数的摘要
%   此处显示详细说明
t = length(T);
Uti = zeros(length(M),1);

variance = zeros(length(T),1);
Pri_loss = zeros(length(T),1);
%l=0;

for i = 1:length(T) 
    variance(T(i)) = sum_sigma_d+sigma_g+ sum(cor_var(M,T(i)));
    Pri_loss(T(i)) = exp(-variance(i));
end

for j = 1:length(M)
    Uti(j) = Cal_Uti_Ind( M(j),M, sum_sigma_d,sigma_g,Pri_loss,S,R_d,r_d);
end

Ave_Uti = mean(Uti);

end

