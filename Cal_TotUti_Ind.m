function [ Total_Uti ] = Cal_TotUti_Ind( M,T, sum_sigma_d,sigma_g,cor_var,S,R_d,r_d )
%CAL_TOTUTI_REPORTER 此处显示有关此函数的摘要
%   此处显示详细说明
t = length(T);
Uti = zeros(t,1);

variance = zeros(length(T),1);
Pri_loss = zeros(length(T),1);
%l=0;

for i = 1:length(T) 
    variance(T(i)) = sum_sigma_d+sigma_g+ sum(cor_var(M,T(i)));
    Pri_loss(T(i)) = exp(-variance(i));
end

for j = 1:t
    Uti(j) = Cal_Uti_Ind( T(j),M, sum_sigma_d,sigma_g,Pri_loss,S,R_d,r_d);
end

Total_Uti = sum(Uti);

end

