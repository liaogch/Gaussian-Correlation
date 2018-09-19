%{
size = 6;
M_min = size;
T = 1:size;

mu_w = 0.5;
sigma_w = 0.4; 
lower_w = 0;
upper_w = inf;
p_w = 0.8;

mu_s = 0.2;
sigma_s = 0.5;
ind_sigma = 0.1;
lower_s = 0;
upper_s = inf;
selfS = 1;
p_s = 0.7;

R_d = 5;
r_d = 0.1;
R_g = 10;
r_a = 0.9;
r_g = 1;


[ S ] = Gen_SocRelMat_TruncNor( size,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s);


[ L , ~ ] = Gen_GauCorM_TruncNor( size,mu_w,sigma_w,lower_w,upper_w,p_w);

cor_var  = Get_Cor_Var( size,L );


r_d = 0.05;
[ M_opt, U_g, Num_reporter,sigma_g ] = ReporterSelection( R_g,r_g,T,S,M_min,cor_var,r_d,r_a); 


[~,~,beta_all] = Get_Threshold( T,T,cor_var,S,r_d );

%}

sigma_g = sigma_g - sigma_g;


iteration = 500;

sigma_d = ones(iteration,size)*1;
step_size = 5;

for t =2:iteration
    for j = 1:size
        temp = 0;
        for i = 1:size
            temp = temp+S(j,i)*exp(-sum(sigma_d(t-1,:))-sigma_g-sum(cor_var(:,i)));
        end
        derivative = -r_d+temp;
       sigma_d(t,j) = max(sigma_d(t-1,j) + step_size*derivative,0);
    end
end




%{
sigma_g = 0;

range = 0:0.01:5;
U_2 = zeros(1,length(range));
U_2_deri = zeros(1,length(range));

for i= 1:length(range)
    sigma_d(1,2) = range(i);
    temp1 = 0;
    for j = 1:size
        temp1 = temp1 + S(2,j)*exp(-sum(sigma_d(1,:))-sigma_g-sum(cor_var(:,j)));
    end
    
    U_2(i) = -r_d * sum(sigma_d(1,:))-temp1;
    
    temp2 = 0;
    for j = 1:size
        temp2 = temp2+S(2,j)*exp(-sum(sigma_d(1,:))-sigma_g-sum(cor_var(:,j)));
    end
    U_2_deri(i) = -r_d +temp2;
end
plot(range,U_2)
figure()
plot(range,U_2_deri)
%}