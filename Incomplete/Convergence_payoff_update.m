B = 5;

stepsize = 0.002;

alpha = 0.5;

r_d = 0.05;

N = 3;

%sigma = [0 0.5 1 1.5 2];
sigma = [ 1  3  5];

sigma_g = 0;

iteration = 500000;

gap = 20;

beta = 5:gap:5+gap*(N-1);

P_0 = ones(length(sigma),1)*0.5;

%{
pd = makedist('Uniform','Lower',1,'Upper',10);

beta = random(pd,N,1);
%}

%P_0 = ones(length(sigma),1)/length(sigma);
%P_0 = [0.4 0.4 0.1 0.1]';
%P_0 = ones(length(sigma),N)/length(sigma);
%P_0(:,1) = [0.6 0.2 0.2  0]';
%P_0(:,2) = [0 0.2 0.2 0.6]';
%P_0(:,3) = [0 0.2 0.2 0.6]';

Perception = cell(N,length(sigma),iteration);


for j = 1:N
    Perception{j} = [P_0 ,zeros(length(sigma),iteration)];
end
%}

%{
for j = 1:N
    P{j} = [P_0(:,j) ,zeros(length(sigma),iteration)];
end
%}
U_max = zeros(N,1);
U_min = zeros(N,1);

for j = 1:N
    [U_max(j), U_min(j)]= Find_Uti_max_min( sigma,sigma_g,beta(j),B,r_d,N );
end

sigma_action = zeros(N,iteration);
normalized_uti = zeros(N,iteration);
for i = 1:iteration
    %stepsize = 1/i;
    for j = 1:N
        Probability = exp(alpha*Perception{j}(:,i))/sum(exp(alpha*Perception{j}(:,i)));
        sigma_action(j,i) = Choose_action (sigma,Probability);
    end
    for j = 1:N
        normalized_uti(j,i) = (Utility_Reporter( sum(sigma_action(:,i)),sigma_g,beta(j),B,r_d )-U_min(j))/(U_max(j)-U_min(j));
        P_temp = Perception{j}(:,i);
        P_temp(find(sigma==sigma_action(j,i))) = normalized_uti(j,i);
        Perception{j}(:,i+1) = Perception{j}(:,i) + stepsize*(P_temp-Perception{j}(:,i));
    end
end



