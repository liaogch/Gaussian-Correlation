
%for the data reporters,
%calculating beta threshold for uniform distribution or other distributions
%and BNE
B = 5;

r_d = 0.1;

%a_g = 2.5;

M = 20;

beta_upper = 3;

opd = makedist('Normal',1.5,5);
pd = truncate(opd,0,beta_upper);

%pd = makedist('Beta',4,2);
%pd = makedist('Uniform',0,beta_upper);
%plot(0:0.01:1,pdf(pd,0:0.01:1))

a_g88 = 3;
beta_threshold88 = zeros(length(a_g88),1);
%a_sum_ave_bay = zeros(length(a_g2),1);
%a_sum_ave_complete = zeros(length(a_g2),1);
%u_bay = zeros(length(a_g2),1);
%u_complete = zeros(length(a_g2),1);

%a_exp = zeros(length(a_g2),1);
%Ma_exp = zeros(length(a_g2),1);

for i = 1:length(a_g88)
   %beta_threshold(i) = Cal_betathreshold_uniform( r_d,a_g(i),M,beta_upper );
   beta_threshold88(i)= Cal_betathreshold_general( r_d,a_g88(i),M,beta_upper,pd);
   %fun = @(x) strategy_BNE(x,beta_threshold(i))*pdf(pd,x);
   %a_exp(i) = integral(fun,0,beta_upper,'ArrayValued',true);
   %Ma_exp(i) = M*a_exp(i);
end

plot(a_g88,beta_threshold88,'LineWidth',2);
%plot(a_g,Ma_exp,'LineWidth',2);
%{
for i = 1:length(a_g)
     iteration = 50000;
     
     beta = random(pd,iteration,M);
     beta_max = zeros(iteration,1);
     a_complete = zeros(iteration,1);
      for j = 1:iteration
        beta_max(j) = max(beta(j,:));
        if log(beta_max(j)/r_d) > a_g(i)
            a_complete(j) = log(beta_max(j)/r_d) - a_g(i);
        end
      end
      a_sum_ave_complete(i)  = mean(a_complete);
end

plot(a_g,Ma_exp,'b',a_g,a_sum_ave_complete,'k','LineWidth',2);
 %}   
%{
for i = 1:length(a_g)
    %beta_threshold(i) = Cal_betathreshold_uniform( r_d,a_g(i),M,beta_upper );
   beta_threshold(i)= Cal_betathreshold_general( r_d,a_g(i),M,beta_upper,pd);
    
    iteration = 50000;

    %pd = makedist('Uniform',0,beta_upper);
    
    beta = random(pd,iteration,M);
    a = zeros(iteration,M);
    u = zeros(iteration,M);
    
    a_complete = zeros(iteration,1);
    beta_max = zeros(iteration,1);
    u_com = zeros(iteration,M);
    
    for j = 1:iteration
        
        for k = 1:M
            if beta(j,k)>beta_threshold(i)
                a(j,k) = log(beta(j,k)/beta_threshold(i));
            end
        end
     
        beta_max(j) = max(beta(j,:));
        if log(beta_max(j)/r_d) > a_g(i)
            a_complete(j) = log(beta_max(j)/r_d) - a_g(i);
        end
        
        for k = 1:M
            a_sum = a_g(i)+sum(a(j,:));
            u(j,k) = B - r_d *(a_sum) - beta(j,k) *exp(-a_sum);
            
            a_sum_com = a_g(i)+a_complete(j);
            u_com(j,k) = B - r_d *(a_sum_com) - beta(j,k) *exp(-a_sum_com);
        end
        
        
    end
    
    a_sum = sum(a.');
    a_sum_ave_bay(i) = mean(a_sum);
    a_sum_ave_complete(i)  = mean(a_complete);
    
    u_bay(i) = mean(mean(u));
    u_complete(i) = mean(mean(u_com));
end

figure
plot(a_g,beta_threshold,'LineWidth',2);

figure

plot(a_g,a_sum_ave_bay,'b*-',a_g,a_sum_ave_complete,'ko-','LineWidth',1.5)

figure

plot(a_g,u_bay,'b*-',a_g,u_complete,'ko-','LineWidth',1.5)
%}

