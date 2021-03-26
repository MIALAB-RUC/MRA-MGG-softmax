function [sigma1,sigma2]=UpdateSigma(f,u,L,Q)
[N,M]=size(f);
temp_u=zeros([N,M]);
for m=1:M
    temp_u(:,m)=circshift(u,L(m));
end
temp_cha_fang=(f-temp_u).^2;
temp_fenzi=temp_cha_fang.*Q;
fenzi=sum(temp_fenzi(:));
fenmu1=sum(Q(:));
sigma1=sqrt(fenzi./(fenmu1+1e-50));
temp_fenzi=temp_cha_fang.*(1-Q);
fenzi=sum(temp_fenzi(:));
fenmu2=M*N-fenmu1;
sigma2=sqrt(fenzi./(fenmu2+1e-50));







% temp_u=gallery('circul',u)';
% for i=1:M
%     tempf=repmat(f(i,:),[N,1]);
%     tempq=squeeze(q(i,:,:))';
%     temp_1=(tempf-temp_u-mu1).^2.*tempq;
%     temp1=sum(temp_1,1);
%     temp_2=(tempf-temp_u-mu2).^2.*(1-tempq);
%     temp2=sum(temp_2,1);
% end
% sigma1=sum(sum(temp1.*w))./(sum(sum(w.*sum(q,3))));
% sigma2=sum(sum(temp2.*w))./(sum(sum(w.*sum((1-q),3))));
% end