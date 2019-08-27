%function [ res ] = LyapunovFun_1( P,r_d,beta,sigma,sigma_g,B, NE )
%LYAPUNOVFUN_1 此处显示有关此函数的摘要
%   此处显示详细说明
beta = [5,25];
X1 = 0:0.01:1;
X2 = 0:0.01:1;
dVdt = zeros(length(X1),length(X2));
for X1_index = 1:length(X1);
    for X2_index = 1:length(X2) 
        %P = [0.1 0.9;0.4 0.6];
        P = [X1(X1_index) 1-X1(X1_index);1-X2(X2_index) X2(X2_index)];

        r_d = 0.05;

        sigma = [1 5];

        sigma_g = 0;
        B = 5;
        NE = [1 2];

        [N,K] = size(P); %N:number of users; K:number of actions
        U_max = zeros(N,1);
        U_min = zeros(N,1);

        for j = 1:N
            [U_max(j), U_min(j)]= Find_Uti_max_min( sigma,sigma_g,beta(j),B,r_d,N );
        end

        Exp_u = zeros(N,1);
        for i = 1:N
            u = 0;
            index = K^N;
            action = ones(N,1);
            for j = 1:index
                sum_action = 0;
                P_temp  = 1;
                for r = 1:N
                    sum_action = sum_action + sigma(action(r));
                    P_temp = P_temp * P(r,action(r)); 
                    if mod(j,K^(r-1)) == 0
                        action(r) = action(r)+1;
                        if action(r) > K;
                            action(r) = mod(action(r),K);
                        end
                    end
                end
                    u = u + (Utility_Reporter(sum_action,sigma_g,beta(i),B,r_d )-U_min(i))/(U_max(i)-U_min(i))*P_temp;
            end
            Exp_u(i) = u;
        end

        NE_u = zeros(N,1);
        for i = 1:N
            u = 0;
    
            index = K^N;
            action = ones(N,1);
    
            P_i = P;
    
            P_NE_i = zeros(1,K);
            P_NE_i(NE(i)) = 1;
            P_i(i,:) = P_NE_i;
    
            for j = 1:index
                sum_action = 0;
                P_temp  = 1;
                for r = 1:N
                    sum_action = sum_action + sigma(action(r));
                    P_temp = P_temp * P_i(r,action(r)); 
                    if mod(j,K^(r-1)) == 0
                        action(r) = action(r)+1;
                        if action(r) > K;
                            action(r) = mod(action(r),K);
                        end
                    end
                end
                    u = u + (Utility_Reporter(sum_action,sigma_g,beta(i),B,r_d )-U_min(i))/(U_max(i)-U_min(i))*P_temp;
            end
            NE_u(i) = u;
        end

        res_temp = 0;
        for i =1:N
            res_temp = res_temp + P(i,NE(i))*(Exp_u(i)-NE_u(i)); 
        end
        dVdt(X1_index,X2_index) = res_temp;
    end
end
[X,Y] = meshgrid(X2,X1);

V = zeros(length(X1),length(X2));
for X1_index = 1:length(X1)
    for X2_index = 1:length(X2)
        V(X1_index,X2_index) = 1-X1(X1_index)+  1-X2(X2_index);
    end
end


%contour(X,Y,dVdt,'ShowText','on');
contourf(X,Y,dVdt,[0 0],'--','ShowText','on');
hold on;
contour(X,Y,V,'ShowText','on');