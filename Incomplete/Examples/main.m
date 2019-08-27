%{
beta_H = 3;
beta_L = 0.5;

P_H = 0.2;
P_L = 0.8;
%}
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
    
    