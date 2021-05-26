function [u] = SubSubUpdateu(f,L,Q,sigma1,sigma2)
M=size(f,2);
for m=1:M
    tempf=f(:,m);
    tempf=circshift(tempf,-L(m));  
    f(:,m)=tempf;
end
fenmu1=sigma1.^2.*(1-Q)+sigma2.^2.*Q;
u=sum(fenmu1.*f,2)./(sum(fenmu1,2)+1e-50);
u(u>1)=1;
u(u<0)=0;
end
