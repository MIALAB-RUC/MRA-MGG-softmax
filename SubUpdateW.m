function [w,w_L]=SubUpdateW(f,u,Q,alpha,sigma1,sigma2)
[N,M]=size(f);
T1=emT(f,u,Q,N);
T2=emT(f,u,1-Q,N);
w1=alpha./(sigma1+eps).*exp(T1)+(1-alpha)./(sigma2+eps).*exp(T2);
w1=bsxfun(@times,w1, 1./(sum(w1,1)+eps));
w=zeros(size(w1));
for m=1:M
    tempw=w1(:,m);
    position= find(tempw==max(tempw));
    w(position(1),m)=1;
end
w_L=find(w~=0);
w_L=mod(w_L-1,N);
w_L=reshape(w_L,[1,M]);
end
function T=emT(f,u,Q,N)
f1=f.*Q;
sqnormf=repmat(sum(abs(f1).^2,1),N,1);
fftu=fft(u);
fftf=fft(f1);
fftQ=fft(Q);
Cf=ifft(bsxfun(@times, conj(fftu), fftf));
sqnormu=ifft(bsxfun(@times, conj(fft(u.^2)), fftQ));
T=2*Cf-sqnormf-sqnormu;
end
