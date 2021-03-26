function  u=MRA_MGG_softmax(X_data,alpha,sigma1,sigma2,lam,r)
u=X_data(:,1);
iiter=12;
q=zeros(size(u));
d=u;
p=zeros(size(u));
%%%%%%%%%%
tau=0.248;
% [w,w_L]=SubUpdateW(X_data,u,alpha,sigma1,sigma2);
% Q=SubUpdateQ(X_data,w_L,u,alpha,sigma1,sigma2);
Q=ones(size(X_data));
for iter=1:5
    [w,w_L]=SubUpdateW(X_data,u,Q,alpha,sigma1,sigma2);
    Q=SubUpdateQ(X_data,w_L,u,alpha,sigma1,sigma2);
    temp=u;
    %update u
    for i=1:iiter
        tempq=q;
        q=myproj(q+tau*mygradient(mydiv(q)-(d-p)),lam/r);
        if sqrt(sum((tempq(:)-q(:)).^2))/sqrt(sum((tempq(:)+eps).^2))< 1e-3
            break;
        end
    end
    u=d-p-mydiv(q);
    u(u>1)=1;
    u(u<0)=0;
    error=sqrt(sum((temp(:)-u(:)).^2))/sqrt(sum((temp(:)+eps).^2));
    if error< 1e-3
        break;
    end
   % fprintf('%.4f\n',error);
    d= SubSubUpdated(X_data,w_L,u,p,r,Q,sigma1,sigma2);
    p=p+u-d;
    if mod(iter,3)==0
    %if iter>=3
     [sigma1,sigma2]=UpdateSigma(X_data,u,w_L,Q);
     alpha=SubUpdateAlpha(Q);
%     pause;
    end
end
end

