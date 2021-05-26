function [X,T,shifts,beta_true] = generate_observations_mix(x, M,alpha, sigma1,sigma2)
%M个观测信号
    x = x(:);
    N = length(x);
    X = zeros(N, M);
    rand('state',0);
    shifts = randi(N, M, 1)'- 1;
    beta_true=zeros(N,1);
    for m = 1 : M
        X(:,m ) = circshift(x, shifts(m));
    end
    for j=1:N
        beta_true(j,1)= numel(find(shifts==j-1))/M;
    end 
%下面对这M个观测信号加入混合高斯噪声
T=zeros(M,N);
T=T(:);
temp=round(alpha*M*N);
rand('seed',1234);
t=randperm(N*M);
X=X(:);
randn('state',0);
v1=sigma1*randn(1,temp);
for n=1:temp
    X(t(n))=X(t(n))+v1(n);  
    T(t(n))=1;   
end
randn('state',0);
v2=sigma2*randn(1,M*N-temp);
for n1=temp+1:M*N
    X(t(n1))=X(t(n1))+v2(n1-temp);   
end
X=reshape(X,[N,M]);
T=reshape(T,[N,M]);
% X(X<0)=0;
% X(X>1)=1;

