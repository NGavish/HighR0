function spreadersAllocation=computeSpreadersAllocation(s,Cij,adultAges,Nadult,vaccinesLeftToDistibute,upperBound)

spreadersAllocation=Nadult*0;

[dummy,groupPriority]=sort(sum(Cij,2),'descend');

ix=1;
while vaccinesLeftToDistibute>0 & ix<=numel(groupPriority)
    jx=groupPriority(ix);
    currVacDist=min([s(jx).*Nadult(jx),upperBound(jx).*Nadult(jx),vaccinesLeftToDistibute]);
    spreadersAllocation(jx)=currVacDist./Nadult(jx);
    vaccinesLeftToDistibute=vaccinesLeftToDistibute-currVacDist;
    ix=ix+1;
end
return

