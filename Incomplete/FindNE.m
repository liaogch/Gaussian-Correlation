function [NE,opt_a_H,BR_a_H,opt_a_L,BR_a_L] = FindNE( a_g,beta_H,P_H,beta_L,P_L,r_d,M,range )

a_H = range;
a_L = range;
opt_a_H = zeros(length(a_H),length(a_L));
BR_a_H = zeros(length(a_H),length(a_L));
opt_a_L = zeros(length(a_H),length(a_L));
BR_a_L = zeros(length(a_H),length(a_L));
NE = [];

for i = 1:length(a_H)
    for j = 1:length(a_L)
        %E_a = P_H * a_H(i) + P_L * a_L(j);
        E_expa = P_H * exp(-a_H(i))+P_L * exp(-a_L(j));
        
        opt_a_H(i,j) = (M-1)*log(E_expa)+log(beta_H/r_d/exp(a_g));
        BR_a_H(i,j) = max(opt_a_H(i,j),0);
        
        opt_a_L(i,j) = (M-1)*log(E_expa)+log(beta_L/r_d/exp(a_g));
        BR_a_L(i,j) = max(opt_a_L(i,j),0);
        
        if abs(a_H(i)-BR_a_H(i,j)) <= 0.001 & abs(a_L(j)-BR_a_L(i,j)) <= 0.001
            NE = [ i j a_H(i) a_L(j) BR_a_H(i,j) BR_a_L(i,j)];
            break;
        end
    end
end


end

