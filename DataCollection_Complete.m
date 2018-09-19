function [ U_g_ave,Ind_TotUti_ave,Num_reporter_ave ] = DataCollection_Complete( t,T,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times)
%DATA_COLLECTION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

U_g = zeros(iteration_times,1);
Total_Uti = zeros(iteration_times,1);
Num_reporter = zeros(iteration_times,1);
%sigma_g = zeros(iteration_times,1);

parfor i = 1:iteration_times

    [ L , ~ ] = Gen_GauCorM_TruncNor( t,mu_w,sigma_w,lower_w,upper_w,p_w);
    [ S ] = Gen_SocRelMat_TruncNor( t,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s);
    cor_var  = Get_Cor_Var( t,L );
    [ M_opt, U_g(i), Num_reporter(i),sigma_g ] = ReporterSelection( R_g,r_g,T,S,M_min,cor_var,r_d,r_a); 
    sum_sigma_d = 0;    
    Total_Uti(i) = Cal_TotUti_Ind( M_opt,T, sum_sigma_d,sigma_g,cor_var,S,R_d,r_d );
end


U_g_ave = mean(U_g);
Ind_TotUti_ave = mean(Total_Uti);
Num_reporter_ave = mean(Num_reporter);

end

