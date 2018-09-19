function [ r ] = Get_SocRelMat( T )
%GET_SOCRELMAT 此处显示有关此函数的摘要
%   此处显示详细说明
if T ==3

    r = [1     1      1.5
         1     1      1
         2.5     1.5    1
        ];
  
%{    
    r = [1     1      1
         1     1      1
         1     1     1
        ];
%}
    

end

if T == 4
    
    r = [1     1      2.5  3.5
         1     1      1     1.5
         2.5   1     1      1
         3.5   1.5   1      1
         ];
     
    
    % r = ones(4);
end
