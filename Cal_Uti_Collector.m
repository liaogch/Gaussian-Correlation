function [ U_g ] = Cal_Uti_Collector( R_g,r_g,r_a,sum_sigma_reporter,sigma_g,M)
%CAL_UTI_COLLECTOR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%sum_sigma_reporter = max(threshold-sigma_g,0);

%R_g = R_g + length(M) * 0.01;

R_g = R_g + length(M) * 0.005;

U_g = R_g-r_g*sum_sigma_reporter-r_a*sigma_g;

end

