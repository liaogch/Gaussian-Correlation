function [ Threshold,reporter,Threshold_all ] = Get_Threshold( M,T,cor_var,s,r_d )

%{
D = zeros(t);
for i = 1:t
    D(i,i) = sum(W(i,:))-W(i,i);
end

L = D-W;

%}
%T = [M N];

%{
for i = 1:t
    for j = [1:i-1 i+1:t]
        cor_var(j,i) = 1/Correlation( L,j,i,t); %cor_var(j,i) = 1/var(x_j|x_i)
    end
end
%}
%{
for j=1:m
    for i = m+1:m+n;
        a = Correlation( L,j,i,t);
        cor_var(j,i) = 1/a;
    end
end
%}


%{
for j=1:length(M)
    for i = length(M)+1:length(M)+length(N);
        a = Correlation( L,j,i,T);
        cor_var(j,i) = 1/a;
    end
end
%}

Threshold_all = zeros(1,length(M));
%{
for j = 1:length(M)
    temp = 0;
    for i = 1:length(N)
        temp = temp + r(M(j),N(i))*exp(-sum(cor_var(M,N(i))));
    end
    Threshold_all(j) = log(temp)-log(r_d);
end
%}
for j = 1:length(M)
    temp = 0;
    for i = 1:length(T);
        temp = temp + s(M(j),T(i))*exp(-sum(cor_var(M,T(i))));
    end
    Threshold_all(j) = max(log(temp)-log(r_d),0);
end


[Threshold,reporter_index] = max(Threshold_all);
%Threshold = max(Threshold,0);
reporter = M(reporter_index);

end

