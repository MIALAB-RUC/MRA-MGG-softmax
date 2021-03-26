function q=myproj(p,lambda)
np=abs(p);
temp=max(np./lambda,1)+1e-10;
q=p./temp;
end