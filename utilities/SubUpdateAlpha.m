function alpha=SubUpdateAlpha(Q)
[N,M]=size(Q);
alpha=sum(Q(:))./(M*N+eps);
end

