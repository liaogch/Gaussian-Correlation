function [ beta_threshold ] = Cal_betathreshold_uniform( r_d,a_g,M,beta_max )
%CAL_BETATHRESHOLDA 此处显示有关此函数的摘要
%   此处显示详细说明
%r_d = 0.05;

%a_g = 2.5;

%M = 5;

%p_beta = 

%beta_max = 1;

fun1 = @(x) 1/beta_max/x;

fun2 = @(x) 1/beta_max;

beta_t = 0.001:0.005:beta_max;



result1 = zeros(length(beta_t),1);
result2 = zeros(length(beta_t),1);

beta_threshold = 0;

for i = 1:length(beta_t)
   % if (beta_t(i)<beta_max)
        result1(i) = integral(fun2,0,beta_t(i),'ArrayValued',true)+beta_t(i)*integral(fun1,beta_t(i),beta_max,'ArrayValued',true);
    %else
     %   result1(i) = 1;
    %end
    result2(i) = power(r_d*exp(a_g)/beta_t(i),1/(M-1));
    if (abs(result1(i)-result2(i)) < 0.005)
        beta_threshold = beta_t(i);
        %break;
    end
end

if beta_threshold == 0;
    beta_threshold = beta_max;
end

%plot(beta_t,result1,beta_t,result2);
%beta_threshold

end

