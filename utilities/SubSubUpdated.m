function  d= SubSubUpdated(f,L,u,p,r,Q,sigma1,sigma2)
M=size(f,2);
for m=1:M
    tempf=f(:,m);
    tempf=circshift(tempf,-L(m));  
    f(:,m)=tempf;
end
fenmu1=sigma1.^2.*(1-Q)+sigma2.^2.*Q;
temp=r.*sigma1.^2+sigma2.^2;
d=(sum(fenmu1.*f,2)+temp.*(u+p))./(sum(fenmu1,2)+temp+1e-50);
d(d>1)=1;
d(d<0)=0;
end