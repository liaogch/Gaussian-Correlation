function [ a_gk,A_t,A,Caseindicator,final_a_g_opt,diff,sum_incomplete ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P )
%for the data collector
% boundaries of the intervals for a_g optimal utility

%{
K = 3;

B = 5;


r_d = 0.1;

r_g = 1;
r_c = 0.9;

M = 20;

%a_g = 0.4;


beta = [6 9 12];



P = [0.5 0.4 0.1];
%}


a_gk = zeros(K,1);
A_t = zeros(K,1);
A = zeros(K,1);
a_g_opt = zeros(K-1,1);
Caseindicator = zeros(K-1,1);
final_a_g_opt = 0;
B1 = zeros(K-1,1);
B2 = zeros(K-1,1);
%diff = 0;

for k = 1:K-1
    temp = 0;
    for i = k+1:K
        temp = temp + P(i)*beta(k)/beta(i);
    end
    a_gk(k) = log(beta(k)*power(sum(P(1:k))+temp,M-1))-log(r_d);
end

a_gk(K) =log(beta(K))-log(r_d);



for k = 1:K-1
    temp = 0;
    for j = k+1:K
        temp = temp+ P(j)*beta(k+1)/beta(j);
    end
    B1(k) = r_g*(sum(P(k+1:K)))*(1+(M-1)*sum(P(1:k))/(sum(P(1:k))+M*temp))-r_c;
    B2(k) = r_g*(sum(P(k+1:K)))*(1+(M-1)*sum(P(1:k))/(sum(P(1:k))+M*temp*beta(k)/beta(k+1)))-r_c;
%{   
    temp = 0;
    for j = k+1:K
        temp = temp+ P(j)*beta(k+1)/beta(j);
    end
    A_t(k+1) = log(M*temp/(r_g*(M-1)*sum(P(1:k))*sum(P(k+1:K))/(r_c-r_g*sum(P(k+1:K)))-sum(P(1:k))));
%}
end


for k=1:K-1
   if (B1(k) >= 0)
       A(1:k+1) = 0;
       
       if (k < K-1)
            for i = k+2: K
                A(i) = log(beta(i)/beta(k+1)); 
            end
       end
       
       a_g_opt(k) = a_gk(k+1);
       final_a_g_opt = a_g_opt(k);
       
       Caseindicator(k) = 1;
       
       
   elseif (B2(k) <= 0)
       
           A(1:k) = 0;
           for i = k+1:K
              A(i) = log(beta(i)/beta(k)); 
           end
           a_g_opt(k) = a_gk(k);
           final_a_g_opt = a_g_opt(k);
           Caseindicator(k) = 2;
           break;
   else
        
        A(1:k) = 0;
        
        temp = 0;
        for j = k+1:K
            temp = temp+ P(j)*beta(k+1)/beta(j);
        end
        
        temp_A_k_1 = log(M*temp/(r_g*(M-1)*sum(P(1:k))*sum(P(k+1:K))/(r_c-r_g*sum(P(k+1:K)))-sum(P(1:k))));
        
        for i = k+1:K
           A(i) =  temp_A_k_1 + log(beta(i)/beta(k+1));
        end
        a_g_opt(k) = log(power(P*exp(-A),M-1)*beta(k+1)/r_d/exp(A(k+1)));
        Caseindicator(k) = 3;
        final_a_g_opt = a_g_opt(k);
        break;
       %{
       temp = sum(P(1:k));
       for j = k+1:K
          temp = temp + P(j)*exp(-A_t(k+1)-log(beta(j)/beta(k+1))); 
       end
       a_g_opt(k) =log( power(temp,M-1)*beta(k+1)/r_d/exp(A_t(k+1)));
      %} 
       
   end
end

sum_incomplete = M*P*A + final_a_g_opt;
sum_complete = log(beta(K)/r_d);
diff = sum_complete - sum_incomplete;

%{
for k = 1:K-1
    if (Caseindicator(k)==3 || Caseindicator(k) == 2)
        final_a_g_opt = a_g_opt(k);
        break;
    else
        final_a_g_opt = a_g_opt(k);
    end  
end
%}


%{
a_g = 0:0.01:6;
NE = zeros(K,length(a_g));
sum_NE = zeros(1,length(a_g));
U = zeros(1,length(a_g));
sum_variance = zeros(1,length(a_g));

a_complete = zeros(1,length(a_g));
U_complete = zeros(1,length(a_g));

sum_variance_complete = zeros(1,length(a_g));




for ii = 1:length(a_g)

action_temp = zeros(K,1);
k_t = 0;
for i = 1:K
    if a_g(ii) <= A_gk(i)
        k_t = i-1;
        break;
    end
end

%k_t

action = zeros(K,1);


for a = 0:0.0001:2
    action_temp = zeros(K,1);
    if (k_t<K)
        action_temp(1:k_t) = 0;
        for i = k_t+1:length(action_temp)
            action_temp(i) = a+log(beta(i)/beta(k_t+1)); 
        end
    end
    if (abs(r_d*exp(a_g(ii)+a)-power(P_beta*exp(-action_temp),M-1)*beta(k_t+1))<= 0.005)
      action= action_temp;
      break;
    end
end

NE(:,ii) = action;
sum_NE(ii) = P_beta*action*M;



U(ii) = B-r_g*M*P_beta*action-r_c*a_g(ii);

sum_variance(ii)  = M*P_beta*action+a_g(ii);

    
if a_g(ii) < A_gk(K)
    a_complete(ii) = A_gk(K)-a_g(ii);
else
    a_complete(ii) = 0;
end

U_complete(ii) = B - r_g*a_complete(ii) - r_c*a_g(ii);
sum_variance_complete(ii) = a_g(ii)+a_complete(ii);

end

U_max = U(1);
a_g_opt = 0;
for i = 1:length(a_g)
    if U_max < U(i)
        U_max = U(i);
        a_g_opt = a_g(i);
    end
end

A_gk
a_g_opt
plot(a_g,U,'b',a_g,U_complete,'k')
figure()
plot(a_g,sum_variance,'b',a_g,sum_variance_complete,'k')
%}
%plot(a_g,sum_NE,'b',a_g,a_complete,'k')

%{
plot(a_g,U,'-','LineWidth',2)
hold on;
plot(a_g,U_complete,'r:','LineWidth',2,'color',[0.7 0.5 0.8]);
%}
%{
grid on;
legend('Incomplete information','Complete information')

set(gca,'FontSize',14)

xlabel('Strategy of the data collector a_g','FontSize',15)
ylabel('Expected utility of the data collector','FontSize',15)
pppppppppppppp
%plot(a_g,NE(1,:),'r',a_g,NE(2,:),'k',a_g,NE(3,:),'b','Linewidth',2)

iteration = 100;
h = ones(K,iteration)/K;
h(:,1) = [0.5 0.4 0.1]';
s = ones(K,1);
iter = 1;
while(iter < iteration)
for i = 1:K
    U_max = -100;
    opt_action_index = 0;
    for a =1:K
        U = B - r_d*action(a) - r_d*(M-1)*h(:,iter)'*action - beta(i)*exp(-action(a))*exp(-a_g)*power(h(:,iter)'*exp(-action),M-1);
        if U>U_max
            U_max = U;
            opt_action_index = a;
        end
    end
    s(i) = opt_action_index;
end
hs = zeros(K,1);

for t = 1:K
    hs(s(t)) = hs(s(t))+P_beta(t);
end
h(:,iter+1) = iter/(iter+1)*h(:,iter)+1/(iter+1)*hs;
iter = iter+1;
end
%}



end

