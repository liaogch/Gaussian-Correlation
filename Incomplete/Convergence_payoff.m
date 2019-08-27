B = 5;

b = 0.002;

r_d = 0.05;

N = 2;

%sigma = [0 0.5 1 1.5 2];
sigma = [ 1 4 ];

sigma_g = 0;

iteration = 500000;

gap = 20;

beta = 5:gap:5+gap*(N-1);

%{
pd = makedist('Uniform','Lower',1,'Upper',10);

beta = random(pd,N,1);
%}

%P_0 = [ 0.1 0.1 0.1 0.7]';
P_0 = [0.5 0.5;0.5 0.5]';
%P_0(1) = [0.1 0.9]';
%P_0(2) = [0.9 0.1]';
%P_0 = ones(length(sigma),1)/length(sigma);
%P_0 = [0.4 0.4 0.1 0.1]';
%P_0 = ones(length(sigma),N)/length(sigma);
%P_0(:,1) = [0.6 0.2 0.2  0]';
%P_0(:,2) = [0 0.2 0.2 0.6]';
%P_0(:,3) = [0 0.2 0.2 0.6]';

P = cell(N,length(sigma),iteration);


for j = 1:N
    P{j} = [P_0(:,j) ,zeros(length(sigma),iteration)];
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
    for j = 1:N
        sigma_action(j,i) = Choose_action (sigma,P{j}(:,i));
    end
    for j = 1:N
        normalized_uti(j,i) = (Utility_Reporter( sum(sigma_action(:,i)),sigma_g,beta(j),B,r_d )-U_min(j))/(U_max(j)-U_min(j));
        e = zeros(length(sigma),1);
        e(find(sigma==sigma_action(j,i))) = 1;
        P{j}(:,i+1) = P{j}(:,i) + b*normalized_uti(j,i)*(e-P{j}(:,i));
    end
end



