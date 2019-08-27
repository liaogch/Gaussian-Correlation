K = 2;

B = 5;

r_d = 0.1;

r_g = 1;
r_c = 0.9;

M = 20;

%a_g = 0.4;
A_gk = zeros(K,1);

beta = [9 50];

delta = beta(1)/beta(2);

P_H = 0.4;

P = [1-P_H P_H];





%B1 = P(2)*r_g*(1+(M-1)*P(1)/(P(1)+M*P(2)))-r_c

%B2 = P(2)*r_g*(1+(M-1)*P(1)/(P(1)+M*P(2)*beta(1)/beta(2)))-r_c

%case 2: A_L = 0  A_H \neq 0


%A_H = log(M*P(2)/((M-1)*P(1)*P(2)*r_g/(r_c-P(2)*r_g)-P(1)))

%A_t = log(beta(2)/beta(1))

%a_g = log(power(P(1)+P(2)*((M-1)*P(1)*P(2)*r_g/(r_c-P(2)*r_g)-P(1))/M/P(2),M-1)*beta(2)/r_d/(M*P(2)/((M-1)*P(1)*P(2)*r_g/(r_c-P(2)*r_g)-P(1))))

%sum_incomplete = M*P(2)*A_H+a_g

%sum_complete = log(beta(2)/r_d)


P_H = 0.001:0.001:0.99;
BB1 = zeros(1,length(P_H));
BB2 = zeros(1,length(P_H));
A_H = zeros(1,length(P_H));
A1 = zeros(1,length(P_H));
A2 = zeros(1,length(P_H));
a_g = zeros(1,length(P_H));
sum_incomplete = zeros(1,length(P_H));
sum_complete = zeros(1,length(P_H));
diff = zeros(1,length(P_H));
t = ones(1,length(P_H));
t1 = ones(1,length(P_H));
for i = 1:length(P_H)
    P = [1-P_H(i), P_H(i)];
    BB1(i) = P(2)*r_g*(1+(M-1)*P(1)/(P(1)+M*P(2)))-r_c;
    BB2(i) = P(2)*r_g*(1+(M-1)*P(1)/(P(1)+M*P(2)*beta(1)/beta(2)))-r_c;
    if (BB1(i) >= 0)
       a_g(i) = log(beta(2)/r_d);
       A1(i) = 0;
       A2(i)=0;
       sum_incomplete(i) = M*(P(1)*A1(i)+P(2)*A2(i))+a_g(i);
       sum_complete(i) = log(beta(2)/r_d);
       diff(i) = sum_complete(i)-sum_incomplete(i);
    elseif (BB2(i) <= 0)
        a_g(i) = (M-1)*log(P(1)+P(2)*delta)+log(beta(1)/r_d);
        A1(i) = 0;
        A2(i) = log(1/delta);
        sum_incomplete(i) = M*(P(1)*A1(i)+P(2)*A2(i))+a_g(i);
        sum_complete(i) = log(beta(2)/r_d);
        diff(i) = sum_complete(i)-sum_incomplete(i);
    else
        a_g(i) = log(power(P(1)+P(2)*((M-1)*P(1)*P(2)*r_g/(r_c-P(2)*r_g)-P(1))/M/P(2),M-1)*beta(2)/r_d/(M*P(2)/((M-1)*P(1)*P(2)*r_g/(r_c-P(2)*r_g)-P(1))));
        A1(i) = 0;
        A2(i) = log(M*P(2)/((M-1)*P(1)*P(2)*r_g/(r_c-P(2)*r_g)-P(1)));
        sum_incomplete(i) = M*(P(1)*A1(i)+P(2)*A2(i))+a_g(i);
        sum_complete(i) = log(beta(2)/r_d);
        diff(i) = sum_complete(i)-sum_incomplete(i);
    end
    t(i) = log(power(M*(P(2)*r_c-P(2)*P(2)*r_g)/(M*P(2)*r_g-M*P(2)*P(2)*r_g+P(2)*r_c-r_c),M*P(2)-1)*power(1-1/M,M-1)*power(1-P(2),M-1)*power(r_c/(r_c-P(2)*r_g),M-1));
    A_H(i) = log(M*P(2)/((M-1)*P(1)*P(2)*r_g/(r_c-P(2)*r_g)-P(1)));
    %a_g = log(power(P(1)+P(2)*((M-1)*P(1)*P(2)*r_g/(r_c-P(2)*r_g)-P(1))/M/P(2),M-1)*beta(2)/r_d/(M*P(2)/((M-1)*P(1)*P(2)*r_g/(r_c-P(2)*r_g)-P(1))));
    %sum_incomplete = M*P(2)*A_H+a_g;
    %sum_complete = log(beta(2)/r_d);
    %if(BB1(i)<0 & BB2(i)>0)
        
        %t(i) = log(power(M*(P(2)*r_c-P(2)*P(2)*r_g)/(M*P(2)*r_g-M*P(2)*P(2)*r_g+P(2)*r_c-r_c),M*P(2)-1)*power(1-1/M,M-1)*power(1-P(2),M-1)*power(r_c/(r_c-P(2)*r_g),M-1));
        
        %M*P(2)*(r_g-r_c)-r_c*(1-P(2))
        %power((1-1/M)*(1-P(2))*r_c/(r_c-P(2)*r_g),M-M*P(2))*power(P(2)*r_c*(M-1)/(M*P(2)*r_g-r_c),M*P(2)-1)
        %(1-P(2))*(M*P(2)*r_g-r_c)/M/P(2)/(r_c-P(2)*r_g)
        %beta(1)/beta(2)
        %power((1-1/M)*(1-P(2))*r_c/(r_c-P(2)*r_g),M-M*P(2))
        %power((1-1/M)*(1-P(2))*r_c/(r_c-P(2)*r_g),M-M*P(2))*power(beta(1)/beta(2)/(1-P(2)+P(2)*beta(1)/beta(2)),1-M*P(2))
        %power(beta(1)/beta(2)/(1-P(2)+P(2)*beta(1)/beta(2)),1-M*P(2))
        %t(i) = log(power(M*(P(2)*r_c-P(2)*P(2)*r_g)/(M*P(2)*r_g-M*P(2)*P(2)*r_g+P(2)*r_c-r_c),M*P(2)-1)*power(1-1/M,M-1)*power(1-P(2),M-1)*power(r_c/(r_c-P(2)*r_g),M-1));
        %t(i) = A_H*(M*P(2)-1)+(M-1)*log(P(1)+P(2)*exp(-A_H));
        %t1(i)  = sum_incomplete-sum_complete;
    %end
end
hold on;
plot(P_H,diff)
%plot(P_H,t)
%{
a_h = 0:0.01:log(beta(2)/beta(1));
t = zeros(length(a_h),1);



for i = 1:length(a_h)
    t(i) = a_h(i)*(M*P(2)-1)+(M-1)*log(P(1)+P(2)*exp(-a_h(i)));
end
plot(a_h,t)
%}




%{
for k = 1:K-1
    temp = 0;
    for i = k+1:K
        temp = temp + P(i)*beta(k)/beta(i);
    end
    A_gk(k) = log(beta(k)*power(sum(P(1:k))+temp,M-1))-log(r_d);
end

A_gk(K) =log(beta(K))-log(r_d);
%}




%{
%case 3: A_L = 0; A_H = log(beta_H/beta/L)


sum_incomplete = M*P(2)*log(beta(2)/beta(1))+log(beta(1)*power(P(1)+P(2)*beta(1)/beta(2),M-1))-log(r_d);

sum_complete = log(beta(2)/r_d);
%}

%{
P_H = 0.05:0.05:0.95;
B = zeros(length(P_H),1);
t = zeros(length(P_H),1);
for i = 1:length(P_H)
    B(i) = (P_H(i)*M-P_H(i)^2*M-1+P_H(i))/(M*P_H(i)-M*P_H(i)*P_H(i));
    t(i) = power(1/B(i),M*P_H(i))*B(i)*power(1-P_H(i)+B(i)*P_H(i),M-1);
end
plot(P_H,B);
figure();
plot(P_H,t);
%}

%{
B1 = zeros(1,length(r_g));
B2 = zeros(1,length(r_g));

for i = 1:length(r_g)
    
B1(i) = P(2)*r_g(i)*(1+(M-1)*P(1)/(P(1)+M*P(2)))-r_c(i);

B2(i) = P(2)*r_g(i)*(1+(M-1)*P(1)/(P(1)+M*P(2)*beta(1)/beta(2)))-r_c(i);
%{
if ( B2(i) <0)
    r_g(i)
    break;

end
%}
end
%}
%{
case 3: A_L = 0; A_H = log(beta_H/beta/L)
sum_incomplete < sum_complete
K = 2;
B = 5;
r_d = 0.1;
r_g = 4;
r_c = 3.9;
M = 20;
beta = [6 9];
P = [0.9 0.1];
%}

%{
case 2 A_L = 0; A_H neq 0
sum_incomplete < sum_complete
K = 2;

B = 5;

r_d = 0.1;

r_g = 4;
r_c = 3.9;

M = 20;

%a_g = 0.4;
A_gk = zeros(K,1);

beta = [1 9];

P_H = 0.6;

P = [1-P_H P_H];

%}

