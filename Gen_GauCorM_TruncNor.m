function [ L , W ] = Gen_GauCorM_TruncNor( t,mu,sigma,lower,upper,p)
%GEN_GAUCORM_TRUNC 此处显示有关此函数的摘要
%   此处显示详细说明
% p 非零概率
D  = zeros(t);
W = zeros(t);
nd=makedist('normal','mu',mu,'sigma',sigma);
td=truncate(nd,lower,upper);
%W = random(td,t,t);
%W = zeros(T);
for i = 1:t 
    flag = 0;
    while(flag == 0)
        W(i,i+1:t) = random(td,1,t-i);
        temp = rand(1,t-i) <= p;
        W(i,i+1:t) = W(i,i+1:t).*temp;
        %{
        for j = i:t
            if rand() <= 1-p
                W(i,j) = 0;
            end
        end
        %}
        if (sum(W(i,:)) ~=0)
                flag = 1;
        end
    end
    W(i+1:t,i) = W(i,i+1:t)';
    %W(i,i) = 0;
    D(i,i) = sum(W(i,:));
end
L = D-W;
end

