function [ U ] = Cal_Uti_Ind( j,M, sum_sigma_d,sigma_g,Pri_loss,S,R_d,r_d)
%CAL_UTI 此处显示有关此函数的摘要
%   此处显示详细说明
%variance = zeros(1,length(T));
%l_i = zeros(1,length(T));
%l=0;
%{
for i = 1:length(T) 
    variance(i) = sum_sigma_d+sigma_g+ sum(cor_var(M,T(i)));
    l_i(i) = exp(-variance(i));
    %l_i(k,i) = log(1/variance(k,i)+1);
    %l = l+S(M(j),T(i))*l_i(i);
    l = l+S(j,T(i))*l_i(i);
    %BR(i) = log(r(1,3)*exp(-sigma_2-sigma_g- sum(cor_var(:,3)))/r_d);
end
%}

l = S(j,:)*Pri_loss;

if any(M==j)
    R_d = R_d+length(M) * 0.01;
    R = R_d - r_d*(sum_sigma_d)-r_d*sigma_g;
else
    R = 0;
end
%alpha = 0.5;
%R(k) = R_max - r_d/alpha*sqrt(2*pi*(sigma_var(k)+sum(sigma_d)+sigma_g));
U = R-l;


