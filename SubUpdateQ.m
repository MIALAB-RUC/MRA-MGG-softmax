function Q=SubUpdateQ(f,L,u,alpha,sigma1,sigma2)
[N,M]=size(f);
temp_u=zeros([N,M]);
for m=1:M
    temp_u(:,m)=circshift(u,L(m));
end
temp_cha_fang=(f-temp_u).^2;
p1=alpha.*Gaussian_pdf(temp_cha_fang,sigma1);
p2=(1-alpha).*Gaussian_pdf(temp_cha_fang,sigma2);
Q=p1./(p1+p2+1e-50);
Q(Q>=.5)=1;
Q(Q<.5)=0;
end
function p=Gaussian_pdf(x,sigma)
p=1./(sigma+eps).*exp(-x.^2./(2*sigma^2+1));
end
% w(w>=.5)=1;
% w(w<.5)=0;
% w=w(:);
% [N,M]=size(f);
% ind=find(w==1);
% ind=mod(ind,N);
% Trans_l=ind-1;
% %问题应该是处在ind上
% tempu=zeros(N,M);
% for t=1:numel(Trans_l)
%     tempu(:,t)=circshift(u,Trans_l(t));
% end
% p1=normpdf(f-tempu,mu1,sigma1);
% p2=normpdf(f-tempu,mu2,sigma2);
% Q=alpha*p1./(alpha*p1+(1-alpha)*p2);
% Q(Q>=.5)=1;
% Q(Q<.5)=0;
