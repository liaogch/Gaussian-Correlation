function [ action ] = Choose_action( sigma,P )
%CHOOSE_ACTION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

num = rand(1);
for i = 1:length(sigma);
    if num < sum(P(1:i))
         action = sigma(i);
         break;
    end
end
end

