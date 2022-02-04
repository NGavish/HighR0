function x=randomizePoint(x,N,Nadult,upperBound,vaccinesLeftToDistibute,delta)
x=min(max(x+delta*(rand(size(x))-0.5),0),upperBound);
if max(x)==0 
    x=ones(size(N));
end
x=vaccinesLeftToDistibute/sum(x.*Nadult)*x;

idxVec=upperBound>0;
idx=find(x>upperBound & upperBound>0);

while numel(idx)>0
    vacExtra=(x(idx)-upperBound(idx))'*Nadult(idx);
    x(idx)=max(upperBound(idx)-1e-15,0);
    idxVec(idx)=false;
    n_idx=find(idxVec);
    x(n_idx)=x(n_idx)+(upperBound(n_idx).*Nadult(n_idx)/sum(upperBound(n_idx).*Nadult(n_idx))*vacExtra)./Nadult(n_idx);
    idx=find(x>upperBound & upperBound>0);
end

return