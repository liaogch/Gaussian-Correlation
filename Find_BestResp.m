function [ BR ] = Find_BestResp( sigma_var,N,j,sigma_d,sigma_g,cor_var,r,R_max,r_d )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明


U = zeros(1,length(sigma_var));
l = zeros(1,length(sigma_var));
R = zeros(1,length(sigma_var));

variance = zeros(length(sigma_var),length(N));
l_i = zeros(length(sigma_var),length(N));

for k=1:length(sigma_var)
    for i = 1:length(N) 
        variance(k,i) = sigma_var(k)+sum(sigma_d)+sigma_g+ sum(cor_var(:,N(i)));
        l_i(k,i) = exp(-variance(k,i));
        %l_i(k,i) = log(1/variance(k,i)+1);
        l(k) = l(k)+r(j,N(i))*l_i(k,i);
        %BR(i) = log(r(1,3)*exp(-sigma_2-sigma_g- sum(cor_var(:,3)))/r_d);
    end
    R(k) = R_max - r_d*(sigma_var(k)+sum(sigma_d))-r_d*sigma_g;
    %alpha = 0.5;
    %R(k) = R_max - r_d/alpha*sqrt(2*pi*(sigma_var(k)+sum(sigma_d)+sigma_g));
    U(k) = R(k)-l(k);
end

%plot(sigma_var,U);
[~, pos] = max(U);
BR = sigma_var(pos);
end

