%for the data reporter 1

T = 3;



W =[
    0     1     0.1
    1     0     2
    0.1   2     0
   ];

r = Get_SocRelMat(T);

weight_1_3 = 0.1;
r(1,3)=1.5;

W(1,3) = weight_1_3;
W(3,1) = weight_1_3;


M = 1:2;
N = 3;
cor_var = zeros(T);



sigma_g = 0:0.01:7;
sigma_1_3 = zeros(1,length(sigma_g));
Threshold = zeros(1,length(sigma_g));
index = zeros(1,length(sigma_g));
for k = 1:length(sigma_g)
    [Threshold(k),index(k)] = Get_Threshold( M,N,T,W,r,r_d );
    if index(k) ==1
        sigma_1_3(k) = max(0,Threshold(k) - sigma_g(k));
    else
        sigma_1_3(k) = 0;
    end
end
hold on;
axis([0 7 -1 7])
plot(sigma_g,sigma_1_3,'LineWidth',2);

%set(gca,'FontSize',12);

