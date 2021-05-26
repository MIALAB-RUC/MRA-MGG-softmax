function dp=mydiv(p)
p1x=repmat(0,[size(p,1), size(p,2)]);
p1x(2:end-1,:)=p(2:end-1,:)-p(1:end-2,:);
p1x(1,:)=p(1,:);
p1x(end,:)=-p(end-1,:);
dp=p1x;
end