function [ action ] = Choose_action( sigma,P )
%CHOOSE_ACTION 此处显示有关此函数的摘要
%   此处显示详细说明

num = rand(1);
for i = 1:length(sigma);
    if num < sum(P(1:i))
         action = sigma(i);
         break;
    end
end
end

