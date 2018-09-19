function [ M_opt, U_opt, remain_num,sigma_opt ] = ReporterSelection( R_g,r_g,T,s,M_min,cor_var,r_d,r_a)
%REPORTERSELECTION 此处显示有关此函数的摘要
%   此处显示详细说明
t = length(T);
%cor_var  = Get_Cor_Var( t,L );

M_all = cell(t-M_min+1,1);
U_all = zeros(t-M_min+1,1);
Reporter_rem = zeros(t-M_min+1,1);

i=1;

M_current = T;
%N_current = [];
M_opt = [];
U_opt = -inf;
sigma_opt = 0;

while (length(M_current) >= M_min)
    
    [ Threshold,reporter ] = Get_Threshold( M_current,T,cor_var,s,r_d );
    sigma_g = Threshold;
    sum_sigma_reporter = max(Threshold-sigma_g,0);
    U = Cal_Uti_Collector( R_g,r_g,r_a,sum_sigma_reporter,sigma_g,M_current);
    
    M_all(i) = {M_current};
    U_all(i) = U;
    Reporter_rem(i) = reporter;
    
    if U > U_opt
        U_opt = U;
        M_opt = M_current;
        sigma_opt = sigma_g;
    end

    i = i+1;
    M_current = M_current(M_current~=reporter);
    %N = [N reporter];
end

remain_num = length(M_opt);

