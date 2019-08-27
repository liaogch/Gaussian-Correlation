function [ output_args ] = Convergence( sigma,N )
%CONVERGENCE 此处显示有关此函数的摘要
%   此处显示详细说明
sigma = [1 2 3 4];

sigma_g = 2;

pd = makedist('Uniform','Lower',1,'Upper',10);

beta = random(pd,N);

end

