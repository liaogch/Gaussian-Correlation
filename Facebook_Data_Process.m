function [ S ] = Facebook_Data_Process( filename )
%FACEBOOK_DATA_PROCESS 此处显示有关此函数的摘要
%   此处显示详细说明
name = [filename '.txt'];
data = importdata(name);
len1 = length(unique(data(:,1)));
len2 = length(unique(data(:,2)));
if len1 ~= len2
    S = 0;
    return;
end
lower = min(data(:,2));
upper = max(data(:,2));

S_temp = zeros(upper - lower);

data_temp = data - lower+1;

for i = 1:length(data_temp(:,1))
    S_temp(data_temp(i,1),data_temp(i,2)) = 1;
    S_temp(data_temp(i,2),data_temp(i,1)) = 1;
end

S_temp(all(S_temp ==0,2),:)=[];
S_temp(:,all(S_temp ==0,1))=[];

S = S_temp;

for i = 1:length(S)
    S(i,i)=1;
end
end

