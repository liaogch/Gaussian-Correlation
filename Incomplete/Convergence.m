function [ output_args ] = Convergence( sigma,N )
%CONVERGENCE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
sigma = [1 2 3 4];

sigma_g = 2;

pd = makedist('Uniform','Lower',1,'Upper',10);

beta = random(pd,N);

end

