function [ U_opt,a_g_opt, NE] = Opt_Collector( beta,P,r_d,M,r_g,r_c)


len = length(beta);

a_g = zeros(1,len);
U = zeros(1,len);

for i = 1:len
    a_g(i) = log(beta(i)*power(sum(P(1:i))+sum(beta(i)*P(i+1:len)./beta(i+1:len)),M-1))-log(r_d);
end

a = zeros(len,len);
for i = 1:len
    for j = (i+1):len
        a(i,j) = 0+log(beta(j)/beta(i));
    end
end

for i = 1:len
    U(i) = -r_g * M* a(i,:)*P' -r_c*a_g(i);
end

[U_opt,index] = max(U);

a_g_opt = a_g(index);
NE = a(index,:);
    

end

