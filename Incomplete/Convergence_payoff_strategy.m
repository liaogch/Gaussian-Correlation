B = 5;

b = 0.1;

r_d = 0.05;

N = 3;

%sigma = [0 0.5 1 1.5 2];

sigma = [ 1 2 3];

sigma_max = max(sigma);
sigma_min = min(sigma);
state = sigma_min * (N-1):sigma(2)-sigma(1):sigma_max * (N-1);
state_size = length(state);

sigma_g = 0;

iteration = 1000000;

gap = 20;

beta = 5:gap:5+gap*(N-1);

%{
pd = makedist('Uniform','Lower',1,'Upper',10);

beta = random(pd,N,1);
%}

P_0 = ones(length(sigma),state_size)/length(sigma);

%{
P_0 = ones(length(sigma),N)/length(sigma);
P_0(:,1) = [0.6 0.2 0.2  0]';
P_0(:,2) = [0 0.2 0.2 0.6]';
P_0(:,3) = [0 0.2 0.2 0.6]';
%}

P = cell(N,iteration,length(sigma),state_size);


for j = 1:N
    P{j,1} = P_0 ;
end


U_max = zeros(N,1);
U_min = zeros(N,1);

for j = 1:N
    [U_max(j), U_min(j)]= Find_Uti_max_min( sigma,sigma_g,beta(j),B,r_d,N );
end

sigma_action = zeros(N,iteration);
normalized_uti = zeros(N,iteration);
last_state_index = ones(N,iteration);

for i = 1:iteration
    for j = 1:N
        sigma_action(j,i) = Choose_action (sigma,P{j,i}(:,last_state_index(j,i)));
    end
    for j = 1:N
        last_state_index(j,i+1) = find(state == sum(sigma_action(:,i))-sigma_action(j,i));
        normalized_uti(j,i) = (Utility_Reporter( sum(sigma_action(:,i)),sigma_g,beta(j),B,r_d )-U_min(j))/(U_max(j)-U_min(j));
        e = zeros(length(sigma),1);
        e(find(sigma==sigma_action(j,i))) = 1;
        P{j,i+1} = P{j,i};
        P{j,i+1}(:,last_state_index(j,i)) = P{j,i}(:,last_state_index(j,i)) + b*normalized_uti(j,i)*(e-P{j,i}(:,last_state_index(j,i)));
    end
    
end

%{
for i = iteration-5000+1:iteration
    for j = 1:N
        sigma_action(j,i) = Choose_action (sigma,P{j,i}(:,last_state_index(j,i)));
    end
    for j = 1:N
        last_state_index(j,i+1) = find(state == sum(sigma_action(:,i))-sigma_action(j,i));
        normalized_uti(j,i) = (Utility_Reporter( sum(sigma_action(:,i)),sigma_g,beta(j),B,r_d )-U_min(j))/(U_max(j)-U_min(j));
        e = zeros(length(sigma),1);
        e(find(sigma==sigma_action(j,i))) = 1;
        P{j,i+1} = P{j,i};
        P{j,i+1}(:,last_state_index(j,i)) = P{j,i}(:,last_state_index(j,i)) + b*normalized_uti(j,i)*(e-P{j,i}(:,last_state_index(j,i)));
    end
    
end
%}



