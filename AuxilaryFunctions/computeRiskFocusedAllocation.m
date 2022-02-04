function riskFocusedAllocation=computeRiskFocusedAllocation(s,Cij,adultAges,Nadult,vaccinesLeftToDistibute,upperBound)

riskFocusedAllocation=Nadult*0;
[dummy,groupPriority]=sort(sum(Cij,2));

ix=1;
while vaccinesLeftToDistibute>0 & ix<=numel(groupPriority)
    jx=groupPriority(ix);
    currVacDist=min([s(jx).*Nadult(jx),upperBound(jx).*Nadult(jx),vaccinesLeftToDistibute]);
    riskFocusedAllocation(jx)=currVacDist./Nadult(jx);
    vaccinesLeftToDistibute=vaccinesLeftToDistibute-currVacDist;
    ix=ix+1;
end
% warning([num2str(vaccinesLeftToDistibute),' were not allocated']);
return

