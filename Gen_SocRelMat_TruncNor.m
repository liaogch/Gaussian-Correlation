function [ S ] = Gen_SocRelMat_TruncNor( t,mu,sigma,ind_sigma,lower,upper,selfS,p )
%GEN_SOCRELMAT_TRUNCNOR 此处显示有关此函数的摘要
%   此处显示详细说明
% p :非零概率
S = eye(t);
if sigma ==0 
    mu_vector = ones(t,1)*mu;
else
    nd1=makedist('normal','mu',mu,'sigma',sigma);
    td1=truncate(nd1,lower,upper);
    mu_vector = random(td1,t,1);
end

for i = 1:t
    nd2=makedist('normal','mu',mu_vector(i),'sigma',ind_sigma);
    td2=truncate(nd2,lower,upper);
    S(i,:) = random(td2,1,t);
end

temp = zeros(t);

for i = 1:t
    temp(i,i+1:t) = rand(1,t-i) <= p;
end

temp = temp +temp';

S = S.* temp;

for i = 1:t
    S(i,i)=selfS;
end

end

