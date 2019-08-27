
%{
Wmin=0;
Wmax=3;
smin=1;
smax=2;
M=1:40;
N=41:50;
T=50;
r_d = 1;

target_M = 1;

[ L50, W50 ] = Gen_GauCorM( T, Wmin,Wmax);

W40 = W50([1:30 41:50],[1:30 41:50]);
W30 = W50([1:20 41:50],[1:20 41:50]);

s50 = Gen_SocRelMat( T,M,N,smin, smax);
s40 = s50([1:30 41:50],[1:30 41:50]);
s30 = s50([1:20 41:50],[1:20 41:50]);

Wave = 1:0.1:3;
Threshold50 = zeros(1,length(Wave));
Threshold40 = zeros(1,length(Wave));
Threshold30 = zeros(1,length(Wave));
W30_new = zeros(30,30,length(Wave));
cor_var_50 = zeros(50,50,length(Wave));
cor_var_30 = zeros(30,30,length(Wave));
index = zeros(1,length(Wave));

for k = 1:length(Wave)
    [ Threshold50(k),index(k)] = CollectorStrategy_DataCorrelation( Wave(k)*50,target_M,1:40,41:50,50,W50,s50,r_d );
end


for k = 1:length(Wave)
    [ Threshold40(k),index(k)] = CollectorStrategy_DataCorrelation( Wave(k)*40,target_M,1:30,31:40,40,W40,s40,r_d );
end

for k = 1:length(Wave)
    [ Threshold30(k),index(k)] = CollectorStrategy_DataCorrelation( Wave(k)*30,target_M,1:20,21:30,30,W30,s30,r_d );
end



plot(Wave,Threshold50,'k',Wave,Threshold40,'b',Wave,Threshold30,'r','LineWidth',2)
legend('M=40','M=30','M=20')
hold on;
%}
%{
load('Corre_50.mat');%L50,W50,s50
M=1:40;
N=41:50;
T=50;
W=W50;
s=s50;
r_d = 1;
R_g=1;
r_g=0.6;
r_a=0.5;
x = 40:-1:20;
Threshold = zeros(1,length(x));
reporter = zeros(1,length(x));
U_g = zeros(1,length(x));
for k=1:length(x)
    [ Threshold(k),reporter(k),reporter_index,cor_var ] = Get_Threshold( M,N,T,W,s,r_d );
    [ U_g(k) ] = Cal_Uti_Collector( R_g,r_g,r_a,Threshold(k),Threshold(k),x(k));
    M(M==reporter(k))=[];
    T=T-1;
    W(reporter_index,:)=[];
    W(:,reporter_index)=[];
    %s(max_index(k),:)=[];
    %s(:,max_index(k))=[];
end
plot(x,U_g,'-o');
%}

%{
load('Corre_50.mat');%L50,W50,s50
MM=1:40;
N=41:50;
t=50;

s=s50;

Wmin=0;
Wmax=3;


target_M = 1;

[ L50, W50 ] = Gen_GauCorM( t, Wmin,Wmax);
W=W50;
%{
D = zeros(1,t);
for i = 1:t
    D(i) = sum(W(i,:))-W(i,i);
end
%}
r_d = 1;
R_g=1;
r_g=0.6;
r_a=0.5;
removed_reporter = 1:40;
Threshold = zeros(1,length(removed_reporter));
reporter = zeros(1,length(removed_reporter));
U_g = zeros(1,length(removed_reporter));
for k=1:length(removed_reporter)
    M=MM([1:removed_reporter(k)-1 removed_reporter(k)+1:40]);
    %M=MM;
    %W=W50([1:removed_reporter(k)-1 removed_reporter(k)+1:50],[1:removed_reporter(k)-1 removed_reporter(k)+1:50]);
    [ Threshold(k),reporter(k),reporter_index,cor_var ] = Get_Threshold( M,40,N,10,t,W,s,r_d );
    [ U_g(k) ] = Cal_Uti_Collector( R_g,r_g,r_a,Threshold(k),Threshold(k),0);
    %s(max_index(k),:)=[];
    %s(:,max_index(k))=[];
end
plot(removed_reporter,U_g,'-o');


[x1,y1]=max(U_g)

for i=1:40
suma(i)=sum(cor_var(i,:));
end
[x2,y2]=min(suma)
%}

%[x3,y3]=max(D(1:40))

%{
Wmin=0;
Wmax=3;
smin=1;
smax=2;
M=1:20;
N=21:50;
T=50;
r_d = 1;

target_M = 1;

[ L50, W50 ] = Gen_GauCorM( T, Wmin,Wmax);

W40 = W50([1:20 21:40],[1:20 21:40]);
W30 = W50([1:20 21:30],[1:20 21:30]);

s50 = Gen_SocRelMat( T,M,N,smin, smax);
s40 = s50([1:20 21:40],[1:20 21:40]);
s30 = s50([1:20 21:30],[1:20 21:30]);

Wave = 1:0.1:3;
Threshold50 = zeros(1,length(Wave));
Threshold40 = zeros(1,length(Wave));
Threshold30 = zeros(1,length(Wave));
W30_new = zeros(30,30,length(Wave));
cor_var_50 = zeros(50,50,length(Wave));
cor_var_30 = zeros(30,30,length(Wave));
index = zeros(1,length(Wave));

for k = 1:length(Wave)
    [ Threshold50(k),index(k)] = CollectorStrategy_DataCorrelation( Wave(k)*50,target_M,1:20,21:50,50,W50,s50,r_d );
end


for k = 1:length(Wave)
    [ Threshold40(k),index(k)] = CollectorStrategy_DataCorrelation( Wave(k)*40,target_M,1:20,21:40,40,W40,s40,r_d );
end

for k = 1:length(Wave)
    [ Threshold30(k),index(k)] = CollectorStrategy_DataCorrelation( Wave(k)*30,target_M,1:20,21:30,30,W30,s30,r_d );
end



plot(Wave,Threshold50,Wave,Threshold40,Wave,Threshold30)
legend('M=30','M=20','M=10')
hold on;
%}

%{
t = 30;
T = 1:t;
%M = 1:20;
%N = [];
M_min = 20;
mu_w = 1;
sigma_w = 0.4; 
lower_w = 0;
upper_w = inf;
p_w = 0.8;

%mu_s = 0.5;
mu_s = 0.5;
sigma_s = 0.5;
ind_sigma = 0.1;
lower_s = 0;
upper_s = inf;
selfS = 1;
p_s = 0.8;

R_d = 5;
r_d = 0.1;
R_g = 10;
r_a = 0.9;
r_g = 1;
%r_g = 1.5;
%{
mu_w = 3.5:1:5.5;
U_g_w_2 = zeros(length(mu_w),1);
Ind_TotUti_w_2 = zeros(length(mu_w),1);
Num_reporter_w_2 = zeros(length(mu_w),1);
Sigma_g_w_2 = zeros(length(mu_w),1);
%}
%{
p_w = 0.2:0.2:1;
U_g_pw = zeros(length(p_w),1);
Ind_TotUti_pw = zeros(length(p_w),1);
Num_reporter_pw = zeros(length(p_w),1);
%}
%U_g_w_3 = zeros(length(mu_w),1);
%Reporter_TotUti_w_3 = zeros(length(mu_w),1);
%{
sigma_s = 0:0.2:1;
U_g_sigmas_3 = zeros(length(sigma_s),1);
Num_reporter_sigmas_3 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_3 = zeros(length(sigma_s),1);
Sigma_g_sigmas_3 = zeros(length(sigma_s),1);
%}

p_s = 0:0.2:1;
U_g_ps = zeros(length(p_s),1);
Ind_TotUti_ps = zeros(length(p_s),1);
Num_reporter_ps = zeros(length(p_s),1);

%{
U_g_sigmas_2 = zeros(length(sigma_s),1);
Reporter_TotUti_sigmas_2 = zeros(length(sigma_s),1);
U_g_sigmas_3 = zeros(length(sigma_s),1);
Reporter_TotUti_sigmas_3 = zeros(length(sigma_s),1);
%}
%parpool(4);
tic;
for j = 1:length(p_s)
    
iteration_times = 500;
U_g = zeros(iteration_times,1);
Total_Uti = zeros(iteration_times,1);
Num_reporter = zeros(iteration_times,1);
sigma_g = zeros(iteration_times,1);
%{
U_g_2 = zeros(iteration_times,1);
Total_Uti_2 = zeros(iteration_times,1);

U_g_3 = zeros(iteration_times,1);
Total_Uti_3 = zeros(iteration_times,1);
%}
parfor i = 1:iteration_times

    [ L , W ] = Gen_GauCorM_TruncNor( t,mu_w,sigma_w,lower_w,upper_w,p_w);
    [ S ] = Gen_SocRelMat_TruncNor( t,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s(j));
    
    cor_var  = Get_Cor_Var( t,L );
    
    [ M_opt, U_g(i), Num_reporter(i),sigma_g(i) ] = ReporterSelection( R_g,r_g,T,S,M_min,cor_var,r_d,r_a);
    
    
    %threshold = Get_Threshold( M_opt,T,cor_var,s,r_d );
    
    %sum_sigma_d_1 = max(threshold-sigma_g,0);
    sum_sigma_d = 0;
    
    Total_Uti(i) = Cal_TotUti_Ind( M_opt,T, sum_sigma_d,sigma_g(i),cor_var,S,R_d,r_d );
    
    %{
    %Complete information
    cor_var  = Get_Cor_Var( t,L );
    [threshold_1,~] = Get_Threshold( M,N,cor_var,S,r_d );
    sigma_g_1 = threshold_1;
    sum_sigma_d_1 = max(threshold_1-sigma_g_1,0);
    [ U_g_1(i) ] = Cal_Uti_Collector( R_g,r_g,r_a,sum_sigma_d_1,sigma_g_1);
    [ Total_Uti_1(i) ] = Cal_TotUti_Reporter( M,T, sum_sigma_d_1,sigma_g_1,cor_var,S,R_d,r_d );
   %}
    %{
    %not awared of S
    Con_S_2 = eye(t)*selfS;
    [Con_threshold_2,~,~ ] = Get_Threshold( M,N,cor_var,Con_S_2,r_d );
    
    [threshold_2,~ ] = Get_Threshold( M,N,cor_var,S,r_d );
    
    sigma_g_2 = Con_threshold_2;
    sum_sigma_d_2 = max(threshold_2-sigma_g_2,0);
    
    [ U_g_2(i) ] = Cal_Uti_Collector( R_g,r_g,r_a,sum_sigma_d_2,sigma_g_2);
    [ Total_Uti_2(i) ] = Cal_TotUti_Reporter( M,T, sum_sigma_d_2,sigma_g_2,cor_var,S,R_d,r_d );
    
    
    
    %Complete S, partial W
    
    Con_W_3 = ones(t)*mu_w;
    Con_D_3 = zeros(t);
    for x = 1:t
        Con_W_3(x,x) = 0;
        Con_D_3(x,x) =sum(Con_W_3(x,:));
    end
    Con_L_3 = Con_D_3 - Con_W_3;
    
    [ Con_threshold_3,~,~ ] = Get_Threshold( M,N,t,Con_L_3,S,r_d );
    
    [threshold_3,reporter,cor_var ] = Get_Threshold( M,N,t,L,S,r_d );
    sigma_g_3 = Con_threshold_3;
    sum_sigma_d_3 = max(threshold_3-sigma_g_3,0);
    
    [ U_g_3(i) ] = Cal_Uti_Collector( R_g,r_g,r_a,sum_sigma_d_3,sigma_g_3);
    [ Total_Uti_3(i) ] = Cal_TotUti_Reporter( M,T, sum_sigma_d_3,sigma_g_3,cor_var,S,R_d,r_d );
    %}
    
end
%{
U_g_sigmas_1(j) = mean(U_g);
Reporter_TotUti_sigmas_1(j) = mean(Total_Uti);
%}
%{
U_g_sigmas_3(j) = mean(U_g);
Ind_TotUti_sigmas_3(j) = mean(Total_Uti);
Num_reporter_sigmas_3(j) = mean(Num_reporter);
Sigma_g_sigmas_3(j) = mean(sigma_g);
%}
%{
U_g_w_2(j) = mean(U_g);
Ind_TotUti_w_2(j) = mean(Total_Uti);
Num_reporter_w_2(j) = mean(Num_reporter);
Sigma_g_w_2(j) = mean(sigma_g);
%}

U_g_ps(j) = mean(U_g);
Ind_TotUti_ps(j) = mean(Total_Uti);
Num_reporter_ps(j) = mean(Num_reporter);

%{
U_g_sigmas_2(j) = mean(U_g_2);
Reporter_TotUti_sigmas_2(j) = mean(Total_Uti_2);

U_g_sigmas_3(j) = mean(U_g_3);
Reporter_TotUti_sigmas_3(j) = mean(Total_Uti_3);
%}
end
toc
%delete(gcp('nocreate'));
figure(1);
hold on;
%plot(sigma_s,U_g_sigmas_1,'-ro',sigma_s,U_g_sigmas_2,'-k^',sigma_s,U_g_sigmas_3,'-bs','LineWidth',2);
plot(p_s,U_g_ps,'-bs','LineWidth',2);
set(gca,'FontSize',12);
% '-ro','-k^','-bs'
figure(2)
hold on;
%plot(sigma_s,Ind_TotUti_sigmas_1,'-ro',sigma_s,Reporter_TotUti_sigmas_2,'-k^',sigma_s,Reporter_TotUti_sigmas_3,'-bs','LineWidth',2);
plot(p_s,Ind_TotUti_ps,'-bs','LineWidth',2);
set(gca,'FontSize',12);
figure(3)
hold on;
%plot(sigma_s,Num_reporter_sigmas_2,'-k^','LineWidth',2);
plot(p_s,Num_reporter_ps,'-bs','LineWidth',2);
set(gca,'FontSize',12);
%}
%{
figure(4)
hold on;
%plot(sigma_s,Num_reporter_sigmas_2,'-k^','LineWidth',2);
plot(sigma_s,Sigma_g_sigmas_3,'-bs','LineWidth',2);
set(gca,'FontSize',12);
%}
%saveresult([sigma_s' U_g_sigmas_1 ],'E:\liao\MATLAB\Gaussian Correlation\data\','TotUti_Reporter_w_3','.txt')


%load('S_FBdata5.mat');

%{
%t = length(S_FBdata5);
S_FBdata = Facebook_Data_Process( 'data6' );
t = length(S_FBdata);
T = 1:t;
%M = 1:20;
%N = [];


M_min = 50;
mu_w = 1;
sigma_w = 0.4; 
lower_w = 0;
upper_w = inf;
p_w = 0.8;

%mu_s = 0.5;
mu_s = 0.5;
sigma_s = 0;
ind_sigma = 0.1;
lower_s = 0;
upper_s = inf;
selfS = 1;
p_s = 0.5;

R_d = 5;
r_d = 0.1;
R_g = 10;
r_a = 0.9;
r_g = 1;

sigma_s = 0:0.5:2;
U_g_sigmas = zeros(length(sigma_s),1);
Num_reporter_sigmas = zeros(length(sigma_s),1);
Ind_TotUti_sigmas = zeros(length(sigma_s),1);

%parpool(4);
tic
for j = 1:length(sigma_s)
    
iteration_times = 500;
U_g = zeros(iteration_times,1);
Total_Uti = zeros(iteration_times,1);
Num_reporter = zeros(iteration_times,1);


parfor i = 1:iteration_times

    [ L , W ] = Gen_GauCorM_TruncNor( t,mu_w,sigma_w,lower_w,upper_w,p_w);
    [ S ] = Gen_SocRelMat_TruncNor( t,mu_s,sigma_s(j),ind_sigma,lower_s,upper_s,selfS,1 ) .*S_FBdata;
    cor_var  = Get_Cor_Var( t,L );
    
    [ M_opt, U_g(i), Num_reporter(i),sigma_g ] = ReporterSelection( R_g,r_g,T,S,M_min,cor_var,r_d,r_a);
    
    
    sum_sigma_d = 0;
    
    Total_Uti(i) = Cal_TotUti_Ind( M_opt,T, sum_sigma_d,sigma_g,cor_var,S,R_d,r_d );
    
end

U_g_sigmas(j) = mean(U_g);
Ind_TotUti_sigmas(j) = mean(Total_Uti);
Num_reporter_sigmas(j) = mean(Num_reporter);

end
toc
%delete(gcp('nocreate'));
figure(1);
hold on;
plot(sigma_s,U_g_sigmas,'-ro','LineWidth',2);
% '-ro','-k^','-bs'
figure(2)
hold on;
plot(sigma_s,Ind_TotUti_sigmas,'-ro','LineWidth',2);

figure(3)
hold on;
plot(sigma_s,Num_reporter_sigmas,'-ro','LineWidth',2);

%saveresult([sigma_s' U_g_sigmas_1 ],'E:\liao\MATLAB\Gaussian Correlation\data\','TotUti_Reporter_w_3','.txt')
%}

%{
t = 30;
T = 1:t;
%M = 1:20;
%N = [];
M_min = 20;
mu_w = 1;
sigma_w = 0.4; 
lower_w = 0;
upper_w = inf;
p_w = 0.8;

%mu_s = 0.5;
mu_s = 0.5;
sigma_s = 0.5;
ind_sigma = 0.1;
lower_s = 0;
upper_s = inf;
selfS = 1;
p_s = 0.8;

R_d = 5;
r_d = 0.1;
R_g = 10;
r_a = 0.9;
r_g = 1;

iteration_times = 1000;

%r_g = 1.5;
%{
mu_w = 3.5:1:5.5;
U_g_w_2 = zeros(length(mu_w),1);
Ind_TotUti_w_2 = zeros(length(mu_w),1);
Num_reporter_w_2 = zeros(length(mu_w),1);
Sigma_g_w_2 = zeros(length(mu_w),1);
%}
%{
p_w = 0.2:0.2:1;
U_g_pw = zeros(length(p_w),1);
Ind_TotUti_pw = zeros(length(p_w),1);
Num_reporter_pw = zeros(length(p_w),1);
%}
%U_g_w_3 = zeros(length(mu_w),1);
%Reporter_TotUti_w_3 = zeros(length(mu_w),1);
%{
sigma_s = 0:0.2:1;
U_g_sigmas_3 = zeros(length(sigma_s),1);
Num_reporter_sigmas_3 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_3 = zeros(length(sigma_s),1);
Sigma_g_sigmas_3 = zeros(length(sigma_s),1);
%}
%{
p_w = 0.4:0.2:1;
U_g_pw_1 = zeros(length(p_w),1);
Ind_TotUti_pw_1 = zeros(length(p_w),1);
Num_reporter_pw_1 = zeros(length(p_w),1);

U_g_pw_2 = zeros(length(p_w),1);
Ind_TotUti_pw_2 = zeros(length(p_w),1);
Num_reporter_pw_2 = zeros(length(p_w),1);

U_g_pw_3 = zeros(length(p_w),1);
Ind_TotUti_pw_3 = zeros(length(p_s),1);
Num_reporter_pw_3 = zeros(length(p_s),1);


U_g_pw_4 = zeros(length(p_w),1);
Ind_TotUti_pw_4 = zeros(length(p_w),1);
Num_reporter_pw_4 = zeros(length(p_w),1);
%}
sigma_s = 0:0.2:1;
U_g_sigmas_1 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_1 = zeros(length(sigma_s),1);
Num_reporter_sigmas_1 = zeros(length(sigma_s),1);

U_g_sigmas_2 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_2 = zeros(length(sigma_s),1);
Num_reporter_sigmas_2 = zeros(length(sigma_s),1);

U_g_sigmas_3 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_3 = zeros(length(sigma_s),1);
Num_reporter_sigmas_3 = zeros(length(sigma_s),1);


U_g_sigmas_4 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_4 = zeros(length(sigma_s),1);
Num_reporter_sigmas_4 = zeros(length(sigma_s),1);

%parpool(4);
tic;
for j = 1:length(sigma_s) 
    [ U_g_sigmas_4(j),Ind_TotUti_sigmas_4(j),Num_reporter_sigmas_4(j) ] = DataCollection_ParS_ComW(t,T,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s(j),ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    [ U_g_sigmas_3(j),Ind_TotUti_sigmas_3(j),Num_reporter_sigmas_3(j) ] = DataCollection_ComS_ParW(t,T,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s(j),ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    [ U_g_sigmas_2(j),Ind_TotUti_sigmas_2(j),Num_reporter_sigmas_2(j) ] = DataCollection_NotSocial_ComW( t,T,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s(j),ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    [ U_g_sigmas_1(j),Ind_TotUti_sigmas_1(j),Num_reporter_sigmas_1(j) ] = DataCollection_Complete( t,T,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s(j),ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
end
toc

figure(1);
hold on;
%plot(sigma_s,U_g_sigmas_1,'-ro',sigma_s,U_g_sigmas_2,'-k^',sigma_s,U_g_sigmas_3,'-bs','LineWidth',2);
plot(mu_w,U_g_w_1,'-ro',mu_w,U_g_w_2,'-k^',mu_w,U_g_w_3,'-bs','LineWidth',2);
set(gca,'FontSize',12);
% '-ro','-k^','-bs','md'
figure(2)
hold on;
%plot(sigma_s,Ind_TotUti_sigmas_1,'-ro',sigma_s,Reporter_TotUti_sigmas_2,'-k^',sigma_s,Reporter_TotUti_sigmas_3,'-bs','LineWidth',2);
plot(mu_w,Ind_TotUti_w_1,'-ro',mu_w,Ind_TotUti_w_2,'-k^',mu_w,Ind_TotUti_w_3,'-bs','LineWidth',2);
set(gca,'FontSize',12);
figure(3)
hold on;
%plot(sigma_s,Num_reporter_sigmas_2,'-k^','LineWidth',2);
plot(mu_w,Num_reporter_w_1,'-ro',mu_w,Num_reporter_w_2,'-k^',mu_w,Num_reporter_w_3,'-bs','LineWidth',2);
set(gca,'FontSize',12);
%}
%{
figure(1);
hold on;
%plot(sigma_s,U_g_sigmas_1,'-ro',sigma_s,U_g_sigmas_2,'-k^',sigma_s,U_g_sigmas_3,'-bs','LineWidth',2);
plot(sigma_s,U_g_sigmas_1,'-ro',sigma_s,U_g_sigmas_2,'-k^',sigma_s,U_g_sigmas_3,'-bs',sigma_s,U_g_sigmas_4,'-md','LineWidth',2);
set(gca,'FontSize',12);
% '-ro','-k^','-bs','md'
figure(2)
hold on;
%plot(sigma_s,Ind_TotUti_sigmas_1,'-ro',sigma_s,Reporter_TotUti_sigmas_2,'-k^',sigma_s,Reporter_TotUti_sigmas_3,'-bs','LineWidth',2);
plot(sigma_s,Ind_TotUti_sigmas_1,'-ro',sigma_s,Ind_TotUti_sigmas_2,'-k^',sigma_s,Ind_TotUti_sigmas_3,'-bs',sigma_s,Ind_TotUti_sigmas_4,'-md','LineWidth',2);
set(gca,'FontSize',12);
figure(3)
hold on;
%plot(sigma_s,Num_reporter_sigmas_2,'-k^','LineWidth',2);
plot(sigma_s,Num_reporter_sigmas_1,'-ro',sigma_s,Num_reporter_sigmas_2,'-k^',sigma_s,Num_reporter_sigmas_3,'-bs',sigma_s,Num_reporter_sigmas_4,'-md','LineWidth',2);
set(gca,'FontSize',12);
%}

%{
filename = 'data6';
M_min = 50;
mu_w = 1;
sigma_w = 0.4; 
lower_w = 0;
upper_w = inf;
p_w = 0.8;

%mu_s = 0.5;
mu_s = 0.5;
sigma_s = 0.5;
ind_sigma = 0.1;
lower_s = 0;
upper_s = inf;
selfS = 1;
p_s = 0.5;

R_d = 5;
r_d = 0.1;
R_g = 10;
r_a = 0.9;
r_g = 1;

%parpool(4);
mu_s = 0.3:0.2:0.9;

U_g_s_1 = zeros(length(mu_s),1);
Ind_TotUti_s_1 = zeros(length(mu_s),1);
Num_reporter_s_1 = zeros(length(mu_s),1);
Rep_AveUti_s_1 = zeros(length(mu_s),1);
Sigma_g_s_1 = zeros(length(mu_s),1);

U_g_s_2 = zeros(length(mu_s),1);
Ind_TotUti_s_2 = zeros(length(mu_s),1);
Num_reporter_s_2 = zeros(length(mu_s),1);
Rep_AveUti_s_2 = zeros(length(mu_s),1);
Sigma_g_s_2 = zeros(length(mu_s),1);

%{
U_g_s_3 = zeros(length(mu_s),1);
Ind_TotUti_s_3 = zeros(length(mu_s),1);
Num_reporter_s_3 = zeros(length(mu_s),1);
Rep_AveUti_s_3 = zeros(length(mu_s),1);
%}

U_g_s_4 = zeros(length(mu_s),1);
Ind_TotUti_s_4 = zeros(length(mu_s),1);
Num_reporter_s_4 = zeros(length(mu_s),1);
Rep_AveUti_s_4 = zeros(length(mu_s),1);
Sigma_g_s_4 = zeros(length(mu_s),1);

iteration_times = 500;
%p_w = 0.4:0.2:1;
tic

for j = 1:length(mu_s)
    mu_w = 0.5;
    [ U_g_s_1(j),Ind_TotUti_s_1(j),Num_reporter_s_1(j),Rep_AveUti_s_1(j),Sigma_g_s_1(j) ] = DataCollection_FBData_Complete(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s(j),sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
end


for j = 1:length(mu_s)
     mu_w = 1.5;
    [ U_g_s_2(j),Ind_TotUti_s_2(j),Num_reporter_s_2(j),Rep_AveUti_s_2(j),Sigma_g_s_2(j) ] = DataCollection_FBData_Complete(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s(j),sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
end

%{
for j = 1:length(mu_s)
     mu_w = 2.5;
    [ U_g_s_3(j),Ind_TotUti_s_3(j),Num_reporter_s_3(j),Rep_AveUti_s_3(j) ] = DataCollection_FBData_Complete(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s(j),sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
end
%}

for j = 1:length(mu_s)
     mu_w = 1;
    [ U_g_s_4(j),Ind_TotUti_s_4(j),Num_reporter_s_4(j),Rep_AveUti_s_4(j),Sigma_g_s_4(j) ] = DataCollection_FBData_Complete(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s(j),sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
end
toc
%delete(gcp('nocreate'));
figure(1);

hold on;
%plot(mu_s,U_g_s_4,'-bs','LineWidth',2);
plot(mu_s,U_g_s_1,'-ro',mu_s,U_g_s_4,'-k^',mu_s,U_g_s_2,'-bs','LineWidth',2);
set(gca,'Xtick',mu_s)
set(gca,'FontSize',12);
% '-ro','-k^','-bs'
%{
figure(2)
hold on;
%plot(mu_s,Ind_TotUti_s_4,'-bs','LineWidth',2);
plot(mu_s,Ind_TotUti_s_1,'-ro',mu_s,Ind_TotUti_s_4,'-k^',mu_s,Ind_TotUti_s_2,'-bs','LineWidth',2);
set(gca,'FontSize',12);
set(gca,'Xtick',mu_s)
%}

figure(3)
hold on;
%plot(mu_s,Num_reporter_s_4,'-bs','LineWidth',2);
plot(mu_s,Num_reporter_s_1,'-ro',mu_s,Num_reporter_s_4,'-k^',mu_s,Num_reporter_s_2,'-bs','LineWidth',2);
set(gca,'FontSize',12);
set(gca,'Xtick',mu_s);

figure(4)
hold on;
%plot(mu_s,Rep_AveUti_s_4,'-bs','LineWidth',2);
plot(mu_s,Rep_AveUti_s_1,'-ro',mu_s,Rep_AveUti_s_4,'-k^',mu_s,Rep_AveUti_s_2,'-bs','LineWidth',2);
set(gca,'FontSize',12);
set(gca,'Xtick',mu_s);

figure(5)
hold on;
%plot(mu_s,Rep_AveUti_s_4,'-bs','LineWidth',2);
plot(mu_s,Sigma_g_s_1,'-ro',mu_s,Sigma_g_s_4,'-k^',mu_s,Sigma_g_s_2,'-bs','LineWidth',2);
set(gca,'FontSize',12);
set(gca,'Xtick',mu_s);


%plot(sigma_s,U_g_sigmas_1,'-ro',sigma_s,U_g_sigmas_2,'-k^',sigma_s,U_g_sigmas_3,'-bs','LineWidth',2);

%plot(mu_w,U_g_w_1,'-ro',mu_w,U_g_w_2,'-k^',mu_w,U_g_w_3,'-bs',mu_w,U_g_w_4,'md','LineWidth',2);
%set(gca,'FontSize',12);
%legend('Scenario 1:complete S and W','Scenario 2:no S; complete W','Scenario 3:complete S; partial W','Scenario 4:partial S; patial W');
% '-ro','-k^','-bs','md'
%}

%filename = 'data6';
filename = 'data4';
M_min = 100;


mu_w = 1;
sigma_w = 0.4; 
lower_w = 0;
upper_w = inf;
p_w = 0.8;

mu_s = 0.5;
sigma_s = 0.5;
%sigma_s = 0.1;
ind_sigma = 0.1;
lower_s = 0;
upper_s = inf;
selfS = 1;
p_s = 0.5;

R_d = 5;
r_d = 0.1;
R_g = 10;
r_a = 0.9;
r_g = 1;



%{
mu_w = 1:0.5:2.5;

U_g_w_1 = zeros(length(mu_w),1);
Ind_TotUti_w_1 = zeros(length(mu_w),1);
Num_reporter_w_1 = zeros(length(mu_w),1);
Rep_AveUti_w_1 = zeros(length(mu_w),1);
Sigma_g_w_1 = zeros(length(mu_w),1);

U_g_w_2 = zeros(length(mu_w),1);
Ind_TotUti_w_2 = zeros(length(mu_w),1);
Num_reporter_w_2 = zeros(length(mu_w),1);
Rep_AveUti_w_2 = zeros(length(mu_w),1);
Sigma_g_w_2 = zeros(length(mu_w),1);
%}
%{
U_g_w_3 = zeros(length(mu_w),1);
Ind_TotUti_w_3 = zeros(length(mu_w),1);
Num_reporter_w_3 = zeros(length(mu_w),1);
Rep_AveUti_w_3 = zeros(length(mu_w),1);
Sigma_g_w_3 = zeros(length(mu_w),1);

U_g_w_4 = zeros(length(mu_w),1);
Ind_TotUti_w_4 = zeros(length(mu_w),1);
Num_reporter_w_4 = zeros(length(mu_w),1);
Rep_AveUti_w_4 = zeros(length(mu_w),1);
Sigma_g_w_4 = zeros(length(mu_w),1);
%}

%sigma_s = 0:0.2:1;
%{
U_g_sigmas_1 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_1 = zeros(length(sigma_s),1);
Num_reporter_sigmas_1 = zeros(length(sigma_s),1);
Rep_AveUti_sigmas_1 = zeros(length(sigma_s),1);
Sigma_g_sigmas_1 = zeros(length(sigma_s),1);


U_g_sigmas_2 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_2 = zeros(length(sigma_s),1);
Num_reporter_sigmas_2 = zeros(length(sigma_s),1);
Rep_AveUti_sigmas_2 = zeros(length(sigma_s),1);
Sigma_g_sigmas_2 = zeros(length(sigma_s),1);
%}
%{
U_g_sigmas_3 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_3 = zeros(length(sigma_s),1);
Num_reporter_sigmas_3 = zeros(length(sigma_s),1);
Rep_AveUti_sigmas_3 = zeros(length(sigma_s),1);
Sigma_g_sigmas_3 = zeros(length(sigma_s),1);



U_g_sigmas_4 = zeros(length(sigma_s),1);
Ind_TotUti_sigmas_4 = zeros(length(sigma_s),1);
Num_reporter_sigmas_4 = zeros(length(sigma_s),1);
Rep_AveUti_sigmas_4 = zeros(length(sigma_s),1);
Sigma_g_sigmas_4 = zeros(length(sigma_s),1);
%}
%p_w = 0.5:0.1:0.9;

%{
U_g_pw_1 = zeros(length(p_w),1);
Ind_TotUti_pw_1 = zeros(length(p_w),1);
Num_reporter_pw_1 = zeros(length(p_w),1);
Rep_AveUti_pw_1 = zeros(length(p_w),1);
Sigma_g_pw_1 = zeros(length(p_w),1);


U_g_pw_2 = zeros(length(p_w),1);
Ind_TotUti_pw_2 = zeros(length(p_w),1);
Num_reporter_pw_2 = zeros(length(p_w),1);
Rep_AveUti_pw_2 = zeros(length(p_w),1);
Sigma_g_pw_2 = zeros(length(p_w),1);
%}


%{
mu_s = 0.3:0.2:0.9;
U_g_s_1 = zeros(length(mu_s),1);
Ind_TotUti_s_1 = zeros(length(mu_s),1);
Num_reporter_s_1 = zeros(length(mu_s),1);
Rep_AveUti_s_1 = zeros(length(mu_s),1);
Sigma_g_s_1 = zeros(length(mu_s),1);


U_g_s_2 = zeros(length(mu_s),1);
Ind_TotUti_s_2 = zeros(length(mu_s),1);
Num_reporter_s_2 = zeros(length(mu_s),1);
Rep_AveUti_s_2 = zeros(length(mu_s),1);
Sigma_g_s_2 = zeros(length(mu_s),1);

U_g_s_3 = zeros(length(mu_s),1);
Ind_TotUti_s_3 = zeros(length(mu_s),1);
Num_reporter_s_3 = zeros(length(mu_s),1);
Rep_AveUti_s_3 = zeros(length(mu_s),1);
Sigma_g_s_3 = zeros(length(mu_s),1);
%}
%{
p_w = 0.5:0.2:0.9;
sigma_s = 0.3:0.2:0.9;
Sigma_g_1 = zeros(length(p_w),length(sigma_s));
%}

mu_s = 0.3:0.2:0.9;
mu_w = 1:0.5:2;

Sigmag_s_w = zeros(length(mu_s),length(mu_w));

iteration_times = 200;


for j = 1:length(mu_s)
    for i = 1:length(mu_w)
    tic
    [ ~,~,~,~,Sigmag_s_w(j,i) ] = DataCollection_FBData_Complete(filename,M_min,mu_w(i),sigma_w,lower_w,upper_w,p_w,mu_s(j),sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ ~,~,~,~,Sigma_g_1(j,i) ] = DataCollection_FBData_Complete(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w(j),mu_s,sigma_s(i),ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_w_1(j),Ind_TotUti_w_1(j),Num_reporter_w_1(j),Rep_AveUti_w_1(j),Sigma_g_w_1(j) ] = DataCollection_FBData_Complete(filename,M_min,mu_w(j),sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_w_2(j),Ind_TotUti_w_2(j),Num_reporter_w_2(j),Rep_AveUti_w_2(j),Sigma_g_w_2(j) ] = DataCollection_FBData_ComS_ParW(filename,M_min,mu_w(j),sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_sigmas_1(j),Ind_TotUti_sigmas_1(j),Num_reporter_sigmas_1(j),Rep_AveUti_sigmas_1(j),Sigma_g_sigmas_1(j) ] = DataCollection_FBData_Complete(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s(j),ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_sigmas_2(j),Ind_TotUti_sigmas_2(j),Num_reporter_sigmas_2(j),Rep_AveUti_sigmas_2(j),Sigma_g_sigmas_2(j) ] = DataCollection_FBData_NotSocial_ComW(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s(j),ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_sigmas_3(j),Ind_TotUti_sigmas_3(j),Num_reporter_sigmas_3(j),Rep_AveUti_sigmas_3(j),Sigma_g_sigmas_3(j) ] = DataCollection_FBData_ComS_ParW(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s(j),ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_sigmas_4(j),Ind_TotUti_sigmas_4(j),Num_reporter_sigmas_4(j),Rep_AveUti_sigmas_4(j),Sigma_g_sigmas_4(j) ] = DataCollection_FBData_ParS_ComW(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s(j),ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_pw_1(j),Ind_TotUti_pw_1(j),Num_reporter_pw_1(j),Rep_AveUti_pw_1(j),Sigma_g_pw_1(j) ] = DataCollection_FBData_Complete(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w(j),mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_pw_2(j),Ind_TotUti_pw_2(j),Num_reporter_pw_2(j),Rep_AveUti_pw_2(j),Sigma_g_pw_2(j) ] = DataCollection_FBData_ComS_ParW(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w(j),mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_s_1(j),Ind_TotUti_s_1(j),Num_reporter_s_1(j),Rep_AveUti_s_1(j),Sigma_g_s_1(j) ] = DataCollection_FBData_Complete(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s(j),sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_s_2(j),Ind_TotUti_s_2(j),Num_reporter_s_2(j),Rep_AveUti_s_2(j),Sigma_g_s_2(j) ] = DataCollection_FBData_NotSocial_ComW(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s(j),sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    %[ U_g_s_3(j),Ind_TotUti_s_3(j),Num_reporter_s_3(j),Rep_AveUti_s_3(j),Sigma_g_s_3(j) ] = DataCollection_FBData_ParS_ComW(filename,M_min,mu_w,sigma_w,lower_w,upper_w,p_w,mu_s(j),sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    toc
    end
end
figure(1)
plot(mu_s,Sigmag_s_w(:,1),'-ro',mu_s,Sigmag_s_w(:,2),'-k^',mu_s,Sigmag_s_w(:,3),'-bs','LineWidth',2);

%{
p_w = 0.8;
sigma_s = 0.5;
sigma_w = 2.5; 

mu_w = 1:0.5:2.5;

U_g_w_11 = zeros(length(mu_w),1);
Ind_TotUti_w_11 = zeros(length(mu_w),1);
Num_reporter_w_11 = zeros(length(mu_w),1);
Rep_AveUti_w_11 = zeros(length(mu_w),1);
Sigma_g_w_11 = zeros(length(mu_w),1);

U_g_w_22 = zeros(length(mu_w),1);
Ind_TotUti_w_22 = zeros(length(mu_w),1);
Num_reporter_w_22 = zeros(length(mu_w),1);
Rep_AveUti_w_22 = zeros(length(mu_w),1);
Sigma_g_w_22 = zeros(length(mu_w),1);

for j = 1:length(mu_w)
    [ U_g_w_11(j),Ind_TotUti_w_11(j),Num_reporter_w_11(j),Rep_AveUti_w_11(j),Sigma_g_w_11(j) ] = DataCollection_FBData_Complete(filename,M_min,mu_w(j),sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
    [ U_g_w_22(j),Ind_TotUti_w_22(j),Num_reporter_w_22(j),Rep_AveUti_w_22(j),Sigma_g_w_22(j) ] = DataCollection_FBData_ComS_ParW(filename,M_min,mu_w(j),sigma_w,lower_w,upper_w,p_w,mu_s,sigma_s,ind_sigma,lower_s,upper_s,selfS,p_s,R_d,r_d,R_g,r_a,r_g,iteration_times);
end



figure(2);
plot(mu_w,U_g_w_11,'-ro',mu_w,U_g_w_22,'-k^','LineWidth',2);
figure(3)
hold on;
plot(mu_w,Ind_TotUti_w_11,'-ro',mu_w,Ind_TotUti_w_22,'-k^','LineWidth',2);
figure(4)
plot(mu_w,Num_reporter_w_11,'-ro',mu_w,Num_reporter_w_22,'-k^','LineWidth',2);
%}

%{
figure(1)
plot(mu_s,U_g_s_1,'-ro','LineWidth',2);
set(gca,'Xtick',mu_s);
figure(2)
plot(mu_s,Ind_TotUti_s_1,'-ro','LineWidth',2);
set(gca,'Xtick',mu_s);
figure(3)
plot(mu_s,Num_reporter_s_1,'-ro','LineWidth',2);
set(gca,'Xtick',mu_s);
figure(4)
plot(mu_s,Sigma_g_s_11,'-ro','LineWidth',2);
set(gca,'Xtick',mu_s);
%}

%{
figure(1)
plot(mu_s,U_g_s_1,'-ro',mu_s,U_g_s_2,'-k^',mu_s,U_g_s_3,'-bs','LineWidth',2);
set(gca,'Xtick',mu_s);
figure(2)
plot(mu_s,Ind_TotUti_s_1,'-ro',mu_s,Ind_TotUti_s_2,'-k^',mu_s,Ind_TotUti_s_3,'-bs','LineWidth',2);
set(gca,'Xtick',mu_s);
figure(3)
plot(mu_s,Num_reporter_s_1,'-ro',mu_s,Num_reporter_s_2,'-k^',mu_s,Num_reporter_s_3,'-bs','LineWidth',2);
set(gca,'Xtick',mu_s);
figure(4)
plot(mu_s,Sigma_g_s_1,'-ro',mu_s,Sigma_g_s_2,'-k^',mu_s,Sigma_g_s_3,'-bs','LineWidth',2);
set(gca,'Xtick',mu_s);
%}

%{
figure(1)
plot(p_w,U_g_pw_1,'-ro',p_w,U_g_pw_2,'-k^','LineWidth',2);
set(gca,'Xtick',p_w);
set(gca,'FontSize',14);
ylabel('Utility of the data collector','FontSize',20);
xlabel('Probability of data correlation p_w','FontSize',20);
figure(2)
plot(p_w,Ind_TotUti_pw_1,'-ro',p_w,Ind_TotUti_pw_2,'-k^','LineWidth',2);
set(gca,'Xtick',p_w);
set(gca,'FontSize',14);
ylabel('Total utility of individuals','FontSize',20);
xlabel('Probability of data correlation p_w','FontSize',20);
figure(3)
plot(p_w,Num_reporter_pw_1,'-ro',p_w,Num_reporter_pw_2,'-k^','LineWidth',2);
set(gca,'Xtick',p_w);
set(gca,'FontSize',14);
ylabel('Size of the reporter set','FontSize',20);
xlabel('Probability of data correlation p_w','FontSize',20);
figure(4)
plot(p_w,Sigma_g_pw_1,'-ro',p_w,Sigma_g_pw_2,'-k^','LineWidth',2);
set(gca,'Xtick',p_w);
set(gca,'FontSize',14);
ylabel('Variance of noise by the data collector','FontSize',20);
xlabel('Probability of data correlation p_w','FontSize',20);

figure(1)
bar(p_w,[U_g_pw_1,U_g_pw_2])
set(gca,'Xtick',p_w);
axis([0.45 0.95 7 9])
set(gca,'FontSize',14);
ylabel('Utility of the data collector','FontSize',20);
xlabel('Data correlation probability p_w','FontSize',20);
legend('Scenario 1: complete S and W','Scenario 3: complete S; partial W');

figure(2)
bar(p_w,[Ind_TotUti_pw_1,Ind_TotUti_pw_2])
axis([0.45 0.95 0 1100])
set(gca,'Xtick',p_w);
axis([0.45 0.95 740 880]);
set(gca,'Ytick',740:20:880);
set(gca,'FontSize',14);
ylabel('Total utility of individuals','FontSize',20);
xlabel('Probability of data correlation p_w','FontSize',20);
legend('Scenario 1: complete S and W','Scenario 3: complete S; partial W');

figure(3)
bar(p_w,[Num_reporter_pw_1,Num_reporter_pw_2])
axis([0.45 0.95 124 142])
set(gca,'Ytick',124:3:142)
set(gca,'FontSize',14);
ylabel('Size of the reporter set','FontSize',20);
xlabel('Data correlation probability p_w','FontSize',20);
legend('Scenario 1: complete S and W','Scenario 3: complete S; partial W');
%}
%{
figure(1)
plot(sigma_s,U_g_sigmas_1,'-ro',sigma_s,U_g_sigmas_2,'-k^',sigma_s,U_g_sigmas_4,'-bs','LineWidth',2);
figure(2)
plot(sigma_s,Ind_TotUti_sigmas_1,'-ro',sigma_s,Ind_TotUti_sigmas_2,'-k^',sigma_s,Ind_TotUti_sigmas_4,'-bs','LineWidth',2);

figure(3)
plot(sigma_s,Num_reporter_sigmas_1,'-ro',sigma_s,Num_reporter_sigmas_2,'-k^',sigma_s,Num_reporter_sigmas_4,'-bs','LineWidth',2);
%}
%{
figure(1);
hold on;
%plot(mu_s,U_g_s_4,'-bs','LineWidth',2);
%plot(mu_w,U_g_w_1,'-ro','LineWidth',2);
plot(mu_w,U_g_w_1,'-ro',mu_w,U_g_w_2,'-k^','LineWidth',2);
set(gca,'Xtick',mu_w)
set(gca,'FontSize',14);
ylabel('Utility of the data collector','FontSize',20);
xlabel('Average data correlation weight \mu_w','FontSize',20);
% '-ro','-k^','-bs','-md'

figure(2)
hold on;
%plot(mu_s,Ind_TotUti_s_4,'-bs','LineWidth',2);
%plot(mu_w,Ind_TotUti_w_1,'-ro','LineWidth',2);
plot(mu_w,Ind_TotUti_w_1,'-ro',mu_w,Ind_TotUti_w_2,'-k^','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',mu_w);
ylabel('Total utility of individuals','FontSize',20);
xlabel('Average data correlation weight \mu_w','FontSize',20);


figure(3)
hold on;

%plot(mu_s,Num_reporter_s_4,'-bs','LineWidth',2);
%plot(mu_w,Num_reporter_w_1,'-ro','LineWidth',2);
plot(mu_w,Num_reporter_w_1,'-ro',mu_w,Num_reporter_w_2,'-k^','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',mu_w);
ylabel('Size of the reporter set','FontSize',20);
xlabel('Average data correlation weight \mu_w','FontSize',20);

%}

%{
figure;
hold on;
plot(sigma_s,U_g_sigmas_1,'-ro','LineWidth',2);
plot(sigma_s,U_g_sigmas_1,'-ro',sigma_s,U_g_sigmas_2,'-k^',sigma_s,U_g_sigmas_4,'-bs','LineWidth',2);
set(gca,'Xtick',sigma_s)
set(gca,'FontSize',14);
ylabel('Utility of the data collector','FontSize',20);
xlabel('Social relationship deviation \sigma_s','FontSize',20);
figure;
plot(sigma_s,Ind_TotUti_sigmas_1,'-ro','LineWidth',2);
plot(sigma_s,Ind_TotUti_sigmas_1,'-ro',sigma_s,Ind_TotUti_sigmas_2,'-k^',sigma_s,Ind_TotUti_sigmas_4,'-bs','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',sigma_s);
ylabel('Total utility of individuals','FontSize',20);
xlabel('Social relationship deviation \sigma_s','FontSize',20);
figure;
plot(sigma_s,Num_reporter_sigmas_1,'-ro','LineWidth',2);
plot(sigma_s,Num_reporter_sigmas_1,'-ro',sigma_s,Num_reporter_sigmas_2,'-k^',sigma_s,Num_reporter_sigmas_4,'-bs','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',sigma_s);
ylabel('Size of the reporter set','FontSize',20);
xlabel('Social relationship deviation \sigma_s','FontSize',20);
%}
%{
%plot(mu_s,U_g_s_4,'-bs','LineWidth',2);
%plot(mu_w,U_g_w_1,'-ro',mu_w,U_g_w_3,'-k^','LineWidth',2);
plot(sigma_s,U_g_sigmas_1,'-ro',sigma_s,U_g_sigmas_2,'-k^',sigma_s,U_g_sigmas_4,'-bs','LineWidth',2);
set(gca,'Xtick',sigma_s)
set(gca,'FontSize',14);
ylabel('Utility of the data collector','FontSize',20);
%xlabel('Average data correlation weight \mu_w','FontSize',20);
xlabel('Social relationship deviation \sigma_s','FontSize',20);
% '-ro','-k^','-bs','-md'

figure(2)
hold on;
%plot(mu_s,Ind_TotUti_s_4,'-bs','LineWidth',2);
plot(sigma_s,Ind_TotUti_sigmas_1,'-ro',sigma_s,Ind_TotUti_sigmas_2,'-k^',sigma_s,Ind_TotUti_sigmas_4,'-bs','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',sigma_s);
ylabel('Total utility of individuals','FontSize',20);
xlabel('Social relationship deviation \sigma_s','FontSize',20);

figure(3)
hold on;
%plot(mu_s,Num_reporter_s_4,'-bs','LineWidth',2);
plot(sigma_s,Num_reporter_sigmas_1,'-ro',sigma_s,Num_reporter_sigmas_2,'-k^',sigma_s,Num_reporter_sigmas_4,'-bs','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',sigma_s);
ylabel('Size of the reporter set','FontSize',20);
xlabel('Social relationship deviation \sigma_s','FontSize',20);

figure(4)
hold on;
%plot(mu_s,Rep_AveUti_s_4,'-bs','LineWidth',2);
plot(sigma_s,Rep_AveUti_sigmas_1,'-ro',sigma_s,Rep_AveUti_sigmas_2,'-k^',sigma_s,Rep_AveUti_sigmas_4,'-bs','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',sigma_s);
ylabel('Average utility of the data reporters','FontSize',20);
xlabel('Social relationship deviation \sigma_s','FontSize',20);

figure(5)
hold on;
%plot(mu_s,Rep_AveUti_s_4,'-bs','LineWidth',2);
plot(sigma_s,Sigma_g_sigmas_1,'-ro',sigma_s,Sigma_g_sigmas_2,'-k^',sigma_s,Sigma_g_sigmas_4,'-bs','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',sigma_s);
ylabel('Variance of noise by the data collector','FontSize',14);
xlabel('Social relationship deviation \sigma_s','FontSize',20);
%}
%{
figure(1);
hold on;
%plot(mu_s,U_g_s_4,'-bs','LineWidth',2);
plot(mu_w(2:6),U_g_w_1(2:6),'-ro',mu_w(2:6),U_g_w_3(2:6),'-k^','LineWidth',2);
set(gca,'Xtick',mu_w(2:6))
set(gca,'FontSize',14);
ylabel('Utility of the data collector','FontSize',20);
xlabel('Average data correlation weight \mu_w','FontSize',20);
% '-ro','-k^','-bs','-md'
figure(2)
hold on;
%plot(mu_s,Ind_TotUti_s_4,'-bs','LineWidth',2);
plot(mu_w(2:6),Ind_TotUti_w_1(2:6),'-ro',mu_w(2:6),Ind_TotUti_w_3(2:6),'-k^','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',mu_w(2:6));
ylabel('Total utility of individuals','FontSize',20);
xlabel('Average data correlation weight \mu_w','FontSize',20);


figure(3)
hold on;

%plot(mu_s,Num_reporter_s_4,'-bs','LineWidth',2);
plot(mu_w(2:6),Num_reporter_w_1(2:6),'-ro',mu_w(2:6),Num_reporter_w_3(2:6),'-k^','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',mu_w(2:6));
ylabel('Size of the reporter set','FontSize',20);
xlabel('Average data correlation weight \mu_w','FontSize',20);
figure(4)
hold on;
%plot(mu_s,Rep_AveUti_s_4,'-bs','LineWidth',2);
plot(mu_w(2:6),Rep_AveUti_w_1(2:6),'-ro',mu_w(2:6),Rep_AveUti_w_3(2:6),'-k^','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',mu_w(2:6));
ylabel('Average utility of the data reporters','FontSize',18);
xlabel('Average data correlation weight \mu_w','FontSize',20);
figure(5)
hold on;
%plot(mu_s,Rep_AveUti_s_4,'-bs','LineWidth',2);
plot(mu_w(2:6),Sigma_g_w_1(2:6),'-ro',mu_w(2:6),Sigma_g_w_3(2:6),'-k^','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'Xtick',mu_w(2:6));
ylabel('Variance of noise by the data collector','FontSize',14);
xlabel('Average data correlation weight \mu_w','FontSize',20);
%}
