function [U_opt,a_g_opt, NE ] = Opt_Collector_p( f,p,r_d,M,r_g,r_c )
pd = makedist('Binomial','N',M-1,'P',p);
P = pdf(pd,0:M-1);
beta = (1:M)*f;
[ U_opt,a_g_opt, NE] = Opt_Collector( beta,P,r_d,M,r_g,r_c);


end

