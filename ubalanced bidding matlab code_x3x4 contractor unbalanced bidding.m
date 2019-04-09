% author: yuhan jiang
% e-mail: yuhan.jiang@marquette.edu

clear
clc
format short;
% (a) Order variables:
M=[0.6834;0.6442]; % Mean
M_ln=[.739,.693]';% Mean of Ln(X)
D_ln=[1.02,.946]'; % Standard Deviations of Ln(X)
%M_ln=[7.39e+005,6.93e+005]';% Mean of Ln(X)
%D_ln=[1.02e+006,9.46e+005]'; % Standard Deviations of Ln(X)
R=[1.0,0.985;0.985,1.0] % corrlection
%R=[1.0,0;0,1.0]
L=chol(R,'lower'); % Lower triangular matrix L=Lo
L_inv = inv(L); % inverse matrix of Lo 

syms x4 x3
g=x4-0.95*x3 % limit state function
gradient_g=gradient(g,[x3,x4]); %  get the gradient of g

% creat array for recording data
x=cell(1,1000); % create array[100] x
z=cell(1,1000); % create array[100] z
u=cell(1,1000); % create array[100] u
grad_g=cell(1,1000); % create array[100] gradient of g
grad_h=cell(1,1000); % create array[100] gradient of h
norm_grad_h=cell(1,1000); % create array[100] norm of the gradient of h(u{})
alpha=cell(1,1000); % create array[100] alpha, direction vector
beta=cell(1,1000); % create array[100] beta, reliability index
h=cell(1,1000); % create array[100] h

%(3a) using the mean value, initial guess for x
x{1}=M;  % initial guess for the x , x1=M
z{1}=[norminv(logncdf(x{1}(1),M_ln(1),D_ln(1)));norminv(logncdf(x{1}(2),M_ln(2),D_ln(2)))]; % x transform to u
u{1}=inv(L)*z{1}; % same as inv(Lo)*z{k} transform to u
k=1;
%repeat ,vpa(x,4) uses at least 4 significant digits
while(k<100)
x3=x{k}(1);% used for calculate gradinet_g
x4=x{k}(2);
grad_g{k}=vpa(subs(gradient_g),4);% evaluated at x{k};

Jux=inv(L)*[lognpdf(x{k}(1),M_ln(1),D_ln(1))/normpdf(z{k}(1)), 0;0, lognpdf(x{k}(2),M_ln(2),D_ln(2))/normpdf(z{k}(2))];
Jxu=inv(Jux);

grad_h{k}= vpa(Jxu'*grad_g{k},4); % equivalent gradient of h(u{k})

norm_grad_h{k}=vpa((grad_h{k}'* grad_h{k})^(1/2),4); % norm of the gradient of h(u{k}) 

alpha{k}=vpa(-1*grad_h{k}/norm_grad_h{k},4);% compute direction vector alpha
beta{k}=vpa(alpha{k}'*u{k},4); % compute reliability index

h{k}=vpa(subs(g),4);% h(u) @k
u{k+1}=vpa(alpha{k}*(beta{k}+h{k}/norm_grad_h{k}),4); % next u
z{k+1}=vpa(L*u{k+1},4);% next z

Pu3=vpa(normcdf(z{k+1}(1),0,1),4);  % P 
Pu4=vpa(normcdf(z{k+1}(2)),4);

xx3=vpa(logninv(Pu3,M_ln(1),D_ln(1)),4); % inverse u to x
xx4=vpa(logninv(Pu4,M_ln(2),D_ln(2)),4);

x{k+1}=[vpa(xx3,4);vpa(xx4,4)];% next x
%(3j) two tolerance value
    if(k>=2 && abs(beta{k}-beta{k-1})<10^(-3) && abs(h{k})<10^(-3))% two tolerance value
         break;
    end
k=k+1;
end
%(4) Present result in table
i=1;
result={'i','x','u','alpha','beta','h'};

while (i<=k)
result=vpa([result;[i,x{i},z{i},alpha{i},beta{i},h{i}]],4);
i=i+1;
end
result%print the result as table.

Sigma=Jxu*Jxu'; % Covariance Matrix
D=diag(sqrt(diag(Sigma)));

numerator = vpa(D*Jux'*alpha{k},4);
denominator = vpa(norm(numerator),4);
gamma = vpa(numerator/denominator,4) % computer gamma
betaFORM=beta{k}
pf1=vpa(normcdf(-1*betaFORM),4)