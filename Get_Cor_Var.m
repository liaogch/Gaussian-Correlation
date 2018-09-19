function [ cor_var ] = Get_Cor_Var( t,L )
%GET_COR_VAR 此处显示有关此函数的摘要
%   此处显示详细说明
cor_var = zeros(t);

for i = 1:t
    for j = i+1:t
        cor_var(j,i) = 1/Correlation( L,j,i,t); %cor_var(j,i) = 1/var(x_j|x_i) symmatric
        %cor_var(i,j) = cor_var(j,i);
    end
end
cor_var = cor_var + cor_var';
end

