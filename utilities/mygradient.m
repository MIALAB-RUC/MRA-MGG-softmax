function gu=mygradient(u)
gu=repmat(0,size(u));
gu(1:end-1,:)=u(2:end,:)-u(1:end-1,:);
end