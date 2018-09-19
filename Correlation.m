function [ a_ji ] = Correlation( L,j,i,T)
%CORRELATION 此处显示有关此函数的摘要
%   此处显示详细说明
%var(x_j|x_i)
l = 1:T;
l(j) = 0;
l(i) = 0;
l = l(l~=0);

W_j = L(j,j);
W_i = L(i,i);
W_l = L(l,l);
W_ji = -L(j,i);
W_jl = -L(j,l);
W_il = -L(i,l);

L_ji = [W_j-W_jl*(inv(W_l))'*W_jl' -W_ji-W_jl*inv(W_l)'*W_il';-W_ji-W_jl*inv(W_l)'*W_il' W_i-W_il*inv(W_l)'*W_il'];
%{
a_ji = W_j;
i = length(l)

for i = 1:length(l)
    a_ji = a_ji - L(j,l(i))^2/L(l(i),l(i));
end
  %}  
a_ji = L_ji(1,1);

end

