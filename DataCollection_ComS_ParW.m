function [ U_g_ave,Ind_TotUti_ave,Num_reporter_ave ] = DataCollection_ComS_ParW(  t,T,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times )
%DATACOLLECTION_COMS_PARW 此处显示有关此函数的摘要
%   此处显示详细说明
U_g = zeros(iteration_times,1);
Total_Uti = zeros(iteration_times,1);
Num_reporter = zeros(iteration_times,1);
%sigma_g = zeros(iteration_times,1);

parfor i = 1:iteration_times
    [ S ] = Gen_SocRelMat_TruncNor( t,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s);
    
    [ L , ~ ] = Gen_GauCorM_TruncNor( t,mu_w,sigma_w,lower_w,upper_w,p_w);
    cor_var  = Get_Cor_Var( t,L );
    
    
    Con_W = ones(t)*mu_w*p_w;
    Con_D = zeros(t);
    for x = 1:t
        Con_W(x,x) = 0;
        Con_D(x,x) =sum(Con_W(x,:));
    end
    Con_L = Con_D - Con_W;
    Con_cor_var = Get_Cor_Var( t,Con_L);
    
    
    
    
    [ Con_M, ~, Num_reporter(i),Con_sigma_g ] = ReporterSelection( R_g,r_g,T,S,M_min,Con_cor_var,r_d,r_a);
    
    [ Threshold,~ ] = Get_Threshold( Con_M,T,cor_var,S,r_d );
    
    %sigma_g_2 = Con_threshold;
    sum_sigma_d = max(Threshold-Con_sigma_g,0);
    
    U_g(i) = Cal_Uti_Collector( R_g,r_g,r_a,sum_sigma_d,Con_sigma_g,Con_M);

    Total_Uti(i) = Cal_TotUti_Ind( Con_M,T, sum_sigma_d,Con_sigma_g,cor_var,S,R_d,r_d );
end


U_g_ave = mean(U_g);
Ind_TotUti_ave = mean(Total_Uti);
Num_reporter_ave = mean(Num_reporter);

end


