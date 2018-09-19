T = 3;
%[ L, W ] = GauCorM( T, Imax);

W = zeros(T);
D = zeros(T);

%mu = [1 2 3 2 1];
%sigma = [1 1.5 2 1.5 1];
%imax = [2 3 4];

W =[
    0     1     0.1
    1     0     3
    0.1   3     0
   ];

weight_1_3 = 0.001:0.001:3;
Threshold = zeros(1,length(weight_1_3));

for k = 1:length(weight_1_3)
W(2,3) = weight_1_3(k);
W(3,2) = weight_1_3(k);


for i = 1:T
    %W(i,:) = mu(i)+sigma(i)^2*randn(1,T);
    %W(i,i+1:T) = Imax(i)*rand(1,T-i);
    %W(i+1:T,i) = W(i,i+1:T)';
    %W(i,i) = 0;
    D(i,i) = sum(W(i,:));
end

L = D-W;

%L = Get_GauCorMat(T);


M = 1:2;
N = 3;
cor_var = zeros(T);

for j=M
    for i = N;
        a = Correlation( L,j,i,T);
        cor_var(j,i) = 1/a;
    end
end

r = Get_SocRelMat(T);
 


R_max = 3;
r_d = 0.001;
%r_d2 = 0.005;
r_g = 0.01;

sigma_g = 1;

Threshold_all = zeros(1,length(M));

for j = 1:length(M)
    temp = 0;
    for i = 1:length(N)
        temp = temp + r(M(j),N(i))*exp(-sum(cor_var(M,N(i))));
    end
    Threshold_all(j) = log(temp)-log(r_d)-sigma_g;
end
Threshold(k) = 0;
index = 0;
[Threshold(k),index] = max(Threshold_all);
Threshold(k) = max(Threshold(k),0);
end

plot(weight_1_3,Threshold);

%{
sigma_2 = 0.2:0.01:2.5;

BR_sigma1 = zeros(1,length(sigma_2));
alpha_1 = zeros(1,length(sigma_2));
for i = 1:length(sigma_2)
    [BR_sigma1(i),alpha_1(i) ]= Cal_BestResp(N,1,[sigma_2(i) sigma_g],cor_var,r,r_d);
end

plot(sigma_2,BR_sigma1);

hold on;

sigma_1 = 0.2:0.01:2.5;
alpha_2 = zeros(1,length(sigma_1));

BR_sigma2 = zeros(1,length(sigma_1));
for i = 1:length(sigma_1)
    [BR_sigma2(i), alpha_2(i)]= Cal_BestResp(N,2,[sigma_1(i) sigma_g],cor_var,r,r_d);
end

plot(BR_sigma2,sigma_1);
%}
%{
sigma_2 = 0:0.1:7;
sigma_1 = 0:0.1:7;
%U = zeros(1,length(sigma_1));
%l_i = zeros(1,length(sigma_1));
%R = zeros(1,length(sigma_1));
BR_sigma1_find  = zeros(1,length(sigma_2));
BR_sigma1_cal  = zeros(1,length(sigma_2));
BR_sigma2_cal  = zeros(1,length(sigma_1));
BR_sigma2_find  = zeros(1,length(sigma_1));
%variance = zeros(1,length(sigma_1));

for i = 1:length(sigma_2)
    BR_sigma1_find(i) = Find_BestResp( 0:0.001:10,N,1,sigma_2(i),sigma_g,cor_var,r,R_max,r_d );
    %BR_sigma1_cal(i) = Cal_BestResp( N,1,[sigma_2(i) sigma_g],cor_var,r,r_d );
end

plot(sigma_2,BR_sigma1_find);


axis([-1 8 -1 8]);
hold on;


for i = 1:length(sigma_1)
    BR_sigma2_find(i) = Find_BestResp( 0:0.001:10,N,2,sigma_1(i), sigma_g,cor_var,r,R_max,r_d );
end

plot(BR_sigma2_find,sigma_1);
xlabel('\sigma^2_2');
ylabel('\sigma^2_1');

Threshold_find = max(BR_sigma1_find(1),BR_sigma2_find(2))



%}

%{
for i=1:length(sigma_1)
    variance(i) = sigma_1(i)+sigma_2+sigma_g+ sum(cor_var(:,3));
    %l_i(i) = exp(-variance(i));
    l_i(i) = log(1/variance(i)+1);
    R(i) = R_max - r_d*(sigma_1(i)+sigma_2)-r_g*sigma_g;
    U(i) = R(i)-r(1,3)*l_i(i);
    %BR(i) = log(r(1,3)*exp(-sigma_2-sigma_g- sum(cor_var(:,3)))/r_d);
end
%}
%BR = Cal_BestResp(N,1,[sigma_2 sigma_g],sum(cor_var(:,N)),r,r_d);

%plot(sigma_1,U);
%}    
        
    
    