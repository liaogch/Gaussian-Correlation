r_c = 1;
r_d = 0.1;
delta = 0.05;
%Scenario 1
s1 = [1 0.9 0.9];
s2 = [0.1 1 0.1];
s3 = [0.1 0.1 1];

c1 = [  1      2];
c2 = [1       2.5];
c3 = [2   2.5   ];


beta_1 = log(s1*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_2 = log(s2*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_3 = log(s3*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_max1 = max(max(beta_1,beta_2),beta_3)
U1 = delta -r_c*beta_max1

c11 = [   1      2];
c22 = [     2.5];
c33 = [  2.5    ];

%beta_11 = log(s1*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_22 = log(s2*[exp(-sum(1./c11)) exp(-sum(1./c22)) exp(-sum(1./c33))]')-log(r_d)
beta_33 = log(s3*[exp(-sum(1./c11)) exp(-sum(1./c22)) exp(-sum(1./c33))]')-log(r_d)
beta_max2 = max(beta_22,beta_33)
U2 = - r_c*beta_max2

select_more_reporter = U1>U2


%Scenario 2

%{
s1 = [1 0.9 0.9];
s2 = [0.1 1 0.1];
s3 = [0.1 0.1 1];
%}

coef = 1.5;
c1 = c1/1.5;
c2 = c2/1.5;
c3 = c3/1.5;


beta_1 = log(s1*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_2 = log(s2*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_3 = log(s3*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_max1 = max(max(beta_1,beta_2),beta_3)
U1 = delta -r_c*beta_max1

c11 = c11/1.5;
c22 = c22/1.5;
c33 = c33/1.5;

%beta_11 = log(s1*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_22 = log(s2*[exp(-sum(1./c11)) exp(-sum(1./c22)) exp(-sum(1./c33))]')-log(r_d)
beta_33 = log(s3*[exp(-sum(1./c11)) exp(-sum(1./c22)) exp(-sum(1./c33))]')-log(r_d)
beta_max2 = max(beta_22,beta_33)
U2 = - r_c*beta_max2

select_more_reporter = U1>U2

betalower = beta_max1 > beta_max2



%Scenario 3
%{
s1 = [1 1 1];
s2 = [0.9 1 0.9];
s3 = [0.9 0.9 1];

c1 = [  1      2];
c2 = [1       2.5];
c3 = [2   2.5   ];


beta_1 = log(s1*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_2 = log(s2*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_3 = log(s3*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_max1 = max(max(beta_1,beta_2),beta_3)
U1 = delta -r_c*beta_max1

c11 = [   1      2];
c22 = [     2.5];
c33 = [  2.5    ];

%beta_11 = log(s1*[exp(-sum(1./c1)) exp(-sum(1./c2)) exp(-sum(1./c3))]')-log(r_d)
beta_22 = log(s2*[exp(-sum(1./c11)) exp(-sum(1./c22)) exp(-sum(1./c33))]')-log(r_d)
beta_33 = log(s3*[exp(-sum(1./c11)) exp(-sum(1./c22)) exp(-sum(1./c33))]')-log(r_d)
beta_max2 = max(beta_22,beta_33)
U2 = - r_c*beta_max2

select_more_reporter = U1>U2
%}


