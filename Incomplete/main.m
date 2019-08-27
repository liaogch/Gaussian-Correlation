%{
beta_H = 3;
beta_L = 0.5;

P_H = 0.2;
P_L = 0.8;
%}
%{
beta_H = 3;
beta_L = 1.5;

P_H = 0.8;
P_L = 0.2;

r_d = 0.1;
M = 5;

a_g_t1 = log(beta_L*power(P_L+P_H*beta_L/beta_H,M-1)/r_d);
a_g_t2 = log (beta_H/r_d);

diff = log (beta_H/beta_L);

range = 0:0.0005:1;


a_g = 0.2:0.1:4;
NE = cell(length(a_g),1);
for i = 1:length(a_g)
    [NE{i,1},opt_a_H,BR_a_H,opt_a_L,BR_a_L] = FindNE( a_g(i),beta_H,P_H,beta_L,P_L,r_d,M,range);
end


r_g = 1;
r_c = 0.9;
U = zeros(length(a_g),1);
a_H = zeros(length(a_g),1);
a_L = zeros(length(a_g),1);
for i =1:length(a_g) 
    if ~isempty(NE{i,1}) 
    a_H(i) = NE{i,1}(1,3);
    a_L(i) = NE{i,1}(1,4);
    U(i) = -r_g*M*[P_H, P_L]*[a_H(i), a_L(i)]'-r_c*a_g(i);
    end
end
plot(a_g,U,'LineWidth',2)
hold on
%}
%{
M = 10;
f = 1;
r_d = 0.1;
r_g = 1;
r_c = 0.9;


beta = zeros(M,1);
P =  zeros(M,1);
p = 0.1:0.005:1;
a_g_opt = zeros(length(p),1);
U_opt = zeros(length(p),1);

NE = cell(length(p),1);



for i = 1:length(p)
    [ U_opt(i),a_g_opt(i), NE{i}] = Opt_Collector_p( f,p(i),r_d,M,r_g,r_c);
end
plot(p,U_opt,'LineWidth',2)

plot(p,a_g_opt,'LineWidth',2)
%}

%{
alpha = 2:0.01:3;
a_g_opt = zeros(length(alpha),1);
U_opt = zeros(length(alpha),1);
NE = cell(length(alpha),1);

beta = zeros(length(alpha),M-1);
P =  zeros(length(alpha),M-1);
for i= 1:length(alpha)
    P(i,:) = [1:M-1].^(-alpha(i));
    P(i,:) = P(i,:)/sum(P(i,:));
    beta(i,:) = [2:M]*f;
    [ U_opt(i),a_g_opt(i), NE{i}] = Opt_Collector( beta(i,:),P(i,:),r_d,M,r_g,r_c);
end

plot(alpha,U_opt,'LineWidth',2);
plot(alpha,a_g_opt,'LineWidth',2)
%}

%{
K = 2;

B = 5;

r_d = 0.1;

r_g = 1;
r_c = 0.9;

M = 10;

%a_g = 0.4;




P2 = 0.001:0.001:0.99;

diff1 = zeros(length(P2),1);
diff2 = zeros(length(P2),1);
diff3 = zeros(length(P2),1);

for i = 1:length(P2)
    
    P = [1-P2(i)  P2(i)];
    beta1 = [6 12];
    [ ~,~,~,~,~,diff1(i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta1,P );
    
    beta2 = [6 18];
    [ ~,~,~,~,~,diff2(i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta2,P );
    
    beta3 = [6 24];
    [ ~,~,~,~,~,diff3(i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta3,P );
end
hold on
plot(P2,diff1,'b-',P2,diff2,'r--',P2,diff3,'k-.','LineWidth',2)
%}

%{
K = 3;

B = 5;

r_d = 0.1;

r_g = 1;
r_c = 0.9;

M = 20;

%a_g = 0.4;


beta = [6 9 12];

P3 = 0.01:0.05:0.96;

P2_axis = 0.01:0.05:0.96;
P3_axis = 0.01:0.05:0.96;


P_cell = cell(length(P2_axis),length(P3));
%diff = zeros(length(P3),1);
diff = ones(length(P2_axis),length(P3_axis))*(-1);
for i = 1:length(P3)
    P2 = 0.01:0.05:1-P3(i)-0.01+0.001;
    for j = 1:length(P2);
        P = [1-P2(j)-P3(i) P2(j)  P3(i)];
        P_cell{j,i} = P;
        [ ~,~,~,~,~,diff(j,i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
    end
end
%{
P2_axis = 0.01:0.05:0.99;
P3_axis = 0.01:0.05:0.99;
%}


s = surf(P3_axis,P2_axis,diff)
set(gca,'zlim',[0,0.65])
view([1 1 1])
xlabel('P_3','FontSize',15)
ylabel('P_2','FontSize',15)
zlabel('S_{com}-S_{in}','FontSize',15)
%}

%{
K = 3;

B = 5;

r_d = 0.1;

r_g = 1;
r_c = 0.9;

M = 3;

%a_g = 0.4;
w_m = 1;

w = w_m + 0.2:0.1:w_m+ 0.8;
sum = zeros(length(w),1);
var_w = zeros(length(w),1);

for i = 1:length(w)
    w_h = w(i);
    
    L_h = [2*w_h -w_h -w_h; -w_h 2*w_h -w_h;-w_h -w_h 2*w_h];
    
    %w_m = 2;
    
    L_m = [2*w_m -w_m -w_m; -w_m 2*w_m -w_m;-w_m -w_m 2*w_m];
    
    w_l = 2*w_m-w(i);
    
    L_l = [2*w_l -w_l -w_l; -w_l 2*w_l -w_l;-w_l -w_l 2*w_l];
    
    beta = [3*exp(-2/Correlation(L_l,1,2,3)) 3*exp(-2/Correlation(L_m,1,2,3)) 3*exp(-2/Correlation(L_h,1,2,3))]

    P = [1/6 2/3 1/6];
    
    var_w(i) = power([w_l w_m w_h] - w_m,2)*P';

    [ a_gk,A_t,A,Caseindicator,final_a_g_opt,diff,sum(i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
end
hold on
%figure();
plot(w,sum,'linewidth',2)
%figure();

%plot(var_w,sum,'linewidth',2)
%}

%{
K = 3;

B = 5;

r_d = 0.1;

r_g = 1;
r_c = 0.9;

M = 3;

%a_g = 0.4;
s_m = 1;

s = s_m+0.1:0.1:s_m+0.8;
sum = zeros(length(s),1);
var_s = zeros(length(s),1);
for i = 1:length(s)

    s_h = s(i);
    
    s_l = 2*s_m-s(i);
    
    beta = [3*s_l*exp(-4/3) 3*s_m*exp(-4/3) 3*s_h*exp(-4/3)]

    P = [1/6 2/3 1/6];
    
    var_s(i) = power([s_l s_m s_h] - s_m,2)*P';
    

    [ a_gk,A_t,A,Caseindicator,final_a_g_opt,diff,sum(i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
end
%figure()
hold on
plot(s,sum,'linewidth',2)
%figure();
hold on
%plot(var_s,sum,'linewidth',2)
%}
%{
K = 3;

B = 5;

r_d = 0.1;

r_g = 1;
r_c = 0.9;

M = 20;

%a_g = 0.4;
w_m = 1;

w =2:0.1:5;
sum = zeros(length(w),1);
var_w = zeros(length(w),1);

for i = 1:length(w)
    w_h = w(i)+0.5;
    
    L_h = [2*w_h -w_h -w_h; -w_h 2*w_h -w_h;-w_h -w_h 2*w_h];
    
    w_m = w(i);
    
    L_m = [2*w_m -w_m -w_m; -w_m 2*w_m -w_m;-w_m -w_m 2*w_m];
    
    w_l = w(i)-0.5;
    
    L_l = [2*w_l -w_l -w_l; -w_l 2*w_l -w_l;-w_l -w_l 2*w_l];
    
    beta = [3*exp(-2/Correlation(L_l,1,2,3)) 3*exp(-2/Correlation(L_m,1,2,3)) 3*exp(-2/Correlation(L_h,1,2,3))];
    
    P = [1/3 1/3 1/3];
    P = [1/4 2/4 1/4];
    P = [1/5 3/5 1/5];
    P = [1/6 4/6 1/6];
    
    
    var_w(i) = power([w_l w_m w_h] - w_m,2)*P';

    [ a_gk,A_t,A,Caseindicator,final_a_g_opt,diff,sum(i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
end
hold on
%figure();
plot(w,sum,':','linewidth',2)
%figure();

%plot(var_w,sum,'linewidth',2)
%}

%{
K = 3;

B = 5;

r_d = 0.1;

r_g = 1;
r_c = 0.9;

M = 20;

%a_g = 0.4;


s = 1:0.1:4;
sum = zeros(length(s),1);
var_s = zeros(length(s),1);
for i = 1:length(s)

    s_h = s(i)+0.5;
    s_m = s(i);
    s_l = s(i)-0.5;
    
    beta = [3*s_l*exp(-4/3) 3*s_m*exp(-4/3) 3*s_h*exp(-4/3)];

    P = [1/3 1/3 1/3];
    P = [1/4 2/4 1/4];
    P = [1/5 3/5 1/5];
    P = [1/6 4/6 1/6];
    
    
    var_s(i) = power([s_l s_m s_h] - s_m,2)*P';
    

    [ a_gk,A_t,A,Caseindicator,final_a_g_opt,diff,sum(i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
end
%figure()
hold on
plot(s,sum,':','linewidth',2)
%figure();
hold on
%plot(var_s,sum,'linewidth',2)

%legend('P=[1/3 1/3 1/3]','P=[1/4 2/4 1/4]','P=[1/5 3/5 1/5]','P=[1/6 4/6 1/6]')
%grid on
%set(gca,'fontsize',15)
%ylabel('S_{in}','fontsize',15)
%xlabel('g_m','fontsize',15)
%}

K = 3;

B = 5;

r_d = 0.1;

r_g = 5;
r_c = 4.9;

M = 8;

beta_ave = 10:5:25;
beta_max = 15:5:30;

diff3 = zeros(length(beta_ave),1);
sum_incomplete3= zeros(length(beta_ave),1);

for i = 1:length(beta_ave)
    d = (beta_max(i)-beta_ave(i))/(K-1)*2;
    
    beta = beta_ave(i)-d:d:beta_ave(i)+d;
    
    P = ones(1,K)/K;
    
    [ ~,~,~,~,~,diff3(i),sum_incomplete3(i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
    
end

K = 5;

diff5 = zeros(length(beta_ave),1);
sum_incomplete5= zeros(length(beta_ave),1);

for i = 1:length(beta_ave)
    d = (beta_max(i)-beta_ave(i))/(K-1)*2;
    
    beta = beta_ave(i)-2*d:d:beta_ave(i)+2*d;
    
    P = ones(1,K)/K;
    
    [ ~,~,~,~,~,diff5(i),sum_incomplete5(i) ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
    
end

K = 7;

diff7 = zeros(length(beta_ave),1);
sum_incomplete7= zeros(length(beta_ave),1);

for i = 1:length(beta_ave)
    d = (beta_max(i)-beta_ave(i))/(K-1)*2;
    
    beta = beta_ave(i)-3*d:d:beta_ave(i)+3*d;
    
    P = ones(1,K)/K;
    
    [ ~,~,~,~,~,diff7(i),sum_incomplete7(i)] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
    
end


plot(beta_ave,sum_incomplete3,'ko-',beta_ave,sum_incomplete5,'r*--',beta_ave,sum_incomplete7,'b^-.','linewidth',2);

%{
hold on
plot(beta_ave,diff3,'ko-',beta_ave,diff5,'r*--',beta_ave,diff7,'b^-.','linewidth',2);
grid on;

legend('K = 3','K = 5','K = 7')
set(gca,'fontsize',14)
ylabel('S_{com}-S_{in}','fontsize',15)
xlabel('Mean of \beta','fontsize',15)
%}
%{
K = 4;

diff4 = zeros(length(beta_ave),1);

for i = 1:length(beta_ave)
    
    d = (beta_max(i)-beta_ave(i))/K*2;
    
    beta = [beta_ave(i)-2*d beta_ave(i)-d beta_ave(i)+d beta_ave(i)+2*d];
    
    P = ones(1,K)/K;
    
    [ ~,~,~,~,~,diff4(i),~ ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
    
end
%}


K = 6;

diff6 = zeros(length(beta_ave),1);

for i = 1:length(beta_ave)
    
    d = (beta_max(i)-beta_ave(i))/K*2;
    
    beta = [beta_ave(i)-3*d beta_ave(i)-2*d beta_ave(i)-d beta_ave(i)+d beta_ave(i)+2*d beta_ave(i)+3*d];
    
    P = ones(1,K)/K;
    
    [ ~,~,~,~,~,diff6(i),~ ] = Ktypes_analytical( K,r_d,r_g,r_c,M,beta,P );
    
end

