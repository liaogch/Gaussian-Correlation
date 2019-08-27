B = 5;

b = 0.01;

r_d = 0.05;

N = 3;

sigma = [0 0.5 1 1.5 2];

sigma_g = 0;

alpha = 0.5;

iteration = 50000;

lambda = 0.1;

c = 1;

%{
beta = 5:4+N;
%}
%{
pd = makedist('Uniform','Lower',1,'Upper',10);

beta = random(pd,N,1);
%}

P_0 = ones(length(sigma),1)/length(sigma);
U_0 = zeros(length(sigma),1)/length(sigma);

P = cell(N,length(sigma),iteration);

U = cell(N,length(sigma),iteration);

for j = 1:N
    P{j} = [P_0 ,zeros(length(sigma),iteration)];
    U{j} = [U_0 ,zeros(length(sigma),iteration)];
end

U_max = zeros(N,1);
U_min = zeros(N,1);

for j = 1:N
    [U_max(j), U_min(j)]= Find_Uti_max_min( sigma,sigma_g,beta(j),B,r_d,N );
end

sigma_action = zeros(N,iteration);
normalized_uti = zeros(N,iteration);
for i = 1:iteration
    lambda = 1/(i+c);
    for j = 1:N
        sigma_action(j,i) = Choose_action (sigma,P{j}(:,i));
    end
    for j = 1:N   
        normalized_uti(j,i) = (Utility_Reporter( sum(sigma_action(:,i)),sigma_g,beta(j),B,r_d )-U_min(j))/(U_max(j)-U_min(j));
        e = zeros(length(sigma),1);
        e(find(sigma==sigma_action(j,i))) = 1;
        for k = 1:length(sigma)
            P{j}(k,i+1) = (1-lambda) * P{j}(k,i) + lambda *exp(alpha*U{j}(k,i))/ sum(exp(alpha*U{j}(:,i)));
            U{j}(k,i+1) = (normalized_uti(j,i)/P{j}(k,i)*e(k)*(normalized_uti(j,i)-U{j}(k,i)) + U{j}(k,i)-U_min(j))/(U_max(j)-U_min(j));
        end
    end
end



