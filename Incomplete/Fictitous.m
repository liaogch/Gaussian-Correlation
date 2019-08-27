K = 3;

B = 2;

r_d = 0.1;

M = 10;

a_g = 0.5;


beta = [2 4 6];

P_beta = [1/3 1/3 1/3];

A_gk = zeros(K,1);

for k = 1:K-1
    temp = 0;
    for i = k+1:K
        temp = temp + P_beta(i)*beta(k)/beta(i);
    end
    A_gk(k) = log(beta(k)*power(sum(P_beta(1:k))+temp,M-1))-log(r_d);
end

A_gk(K) =log(beta(K))-log(r_d);

A_gk

action_temp = zeros(K,1);
k_t = 0;
for i = 1:K
    if a_g <= A_gk(i)
        k_t = i-1;
        break;
    end
end

k_t

action = zeros(K,1);

for a = 0:0.001:2
    action_temp = zeros(K,1);
    if (k_t<K)
        action_temp(1:k_t) = 0;
        for i = k_t+1:length(action_temp)
            action_temp(i) = a+log(beta(i)/beta(k_t+1)); 
        end
    end
    if (abs(r_d*exp(a_g+a)-power(P_beta*exp(-action_temp),M-1)*beta(k_t+1))<= 0.01)
      action= action_temp;
      break;
    end
end
action

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

s