%for the data collector
T = 3;

W =[
    0     1     0.1
    1     0     2
    0.1   2     0
   ];

r = Get_SocRelMat(T);

%weight_1_3 = 0.001:0.001:3;
r_1_3=0:0.01:2;

min = 0;
max = 2;
W=min + (max-min)*rand(20,20); 


M = 1:2;
N = 3;

%{
Threshold = zeros(1,length(weight_1_3));
index = zeros(1,length(weight_1_3));

for k = 1:length(weight_1_3)
    %W(2,3) = weight_1_3(k);
    %W(3,2) = weight_1_3(k);
    [Threshold(k),index(k)] = Get_Threshold( M,N,T,W,r,r_d );
end
hold on;
axis([0 3 -1 6])
plot(weight_1_3,Threshold,'LineWidth',2);
%}

Threshold = zeros(1,length(r_1_3));
index = zeros(1,length(r_1_3));

for k = 1:length(r_1_3)
    r(1,3)=r_1_3(k);
    [Threshold(k),index(k)] = Get_Threshold( M,N,T,W,r,r_d );
end
hold on;
plot(r_1_3,Threshold,'LineWidth',2);

%set(gca,'FontSize',12);

