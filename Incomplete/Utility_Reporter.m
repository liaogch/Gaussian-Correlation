function [ Uti ] = Utility_Reporter( sigma_sum,sigma_g,beta,B,r_d )
%UTILITY_REPORTER �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

Uti = B - r_d * (sigma_sum+sigma_g) - exp(-sigma_sum-sigma_g)*beta;

end

