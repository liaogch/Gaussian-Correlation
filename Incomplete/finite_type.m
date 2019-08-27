%for the data collector
% boundaries of the intervals for a_g optimal utility

K = 3;

B = 5;

r_d = 0.1;

r_g = 1;
r_c = 0.9;

M = 20;

%a_g = 0.4;


beta = [6 9 12];


P = [0.5 0.4 0.1];

A_gk = zeros(K,1);

for k = 1:K-1
    temp = 0;
    for i = k+1:K
        temp = temp + P(i)*beta(k)/beta(i);
    end
    A_gk(k) = log(beta(k)*power(sum(P(1:k))+temp,M-1))-log(r_d);
end

A_gk(K) =log(beta(K))-log(r_d);


a_g = 0:0.01:5;
NE = zeros(K,length(a_g));
sum_NE = zeros(1,length(a_g));
U = zeros(1,length(a_g));

a_complete = zeros(1,length(a_g));
U_complete = zeros(1,length(a_g));

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
    if (abs(r_d*exp(a_g(ii)+a)-power(P*exp(-action_temp),M-1)*beta(k_t+1))<= 0.05)
      action= action_temp;
      break;
    end
end

NE(:,ii) = action;
sum_NE(ii) = P*action*M;



U(ii) = B-r_g*M*P*action-r_c*a_g(ii);

    
if a_g(ii) < A_gk(K)
    a_complete(ii) = A_gk(K)-a_g(ii);
else
    a_complete(ii) = 0;
end

U_complete(ii) = B - r_g*a_complete(ii) - r_c*a_g(ii);

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
plot(a_g,U,'-','LineWidth',2)
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
%}
%plot(a_g,NE(1,:),'r',a_g,NE(2,:),'k',a_g,NE(3,:),'b','Linewidth',2)
%{
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
