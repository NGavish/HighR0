function [vacOfMoreSusceptible,vacOfLessSusceptible,optimalAllocationRes,optimalAllocation,attackRatevacOfMoreSusceptible,attackRatevacOfLessSusceptible]=outcomeAsFunctionofR(Rvalues,Cij,Ni,vaccineCoverage,betaVaccine)

%% Define
defineColors
N=numel(Ni);
populationSize=sum(Ni);
adultAges=1:N;
e=ones(size(Ni));
vaccinesLeftToDistibute=populationSize*vaccineCoverage/100;
spreadersAllocation=computeSpreadersAllocation(e,Cij,adultAges,Ni,populationSize*vaccineCoverage/100,e);
riskFocusedAllocation=computeRiskFocusedAllocation(e,Cij,adultAges,Ni,populationSize*vaccineCoverage/100,e);

%% Optimization problem
problem.options = optimoptions('fmincon','MaxFunctionEvaluations',5e4,'ConstraintTolerance',1e-6,'StepTolerance',1e-10,'Display','none','Algorithm','sqp');%,'Display','iter');
problem.solver = 'fmincon';
problem.Aeq=Ni';problem.Beq=populationSize*vaccineCoverage/100;
problem.A=[];problem.B=[];
problem.lb=0*e;
problem.ub=e;

%% Compute final size as a function of R
M=numel(Rvalues);

optimalAllocationRes=zeros(1,M);
vacOfLessSusceptible=zeros(1,M);
vacOfMoreSusceptible=zeros(1,M);
attackRatevacOfLessSusceptible=zeros(N,M);
attackRatevacOfMoreSusceptible=zeros(N,M);
optimalAllocation=cell(1,M);
[V,d]=eig(Cij);
r=0*e;%0.2*abs(V(:,1))/norm(V(:,1));
for ix=1:M
    R0=Rvalues(ix)
    [vacOfMoreSusceptible(ix),dummy,dummy,data]=computeFinalSize_generalized(spreadersAllocation,adultAges,e,0,R0,Cij,Ni,r,0*e,0*e,0*e,1,betaVaccine,0,1,1);
    attackRatevacOfMoreSusceptible(:,ix)=data(:,2)+data(:,4);
    [vacOfLessSusceptible(ix),dummy,dummy,data]=computeFinalSize_generalized(riskFocusedAllocation,adultAges,e,0,R0,Cij,Ni,r,0*e,0*e,0*e,1,betaVaccine,0,1,1);
    attackRatevacOfLessSusceptible(:,ix)=data(:,2)+data(:,4);

    [optimalAllocation{ix},optimalAllocationRes(ix),exitflag,output]=fmincon(@(x) computeFinalSize_generalized(x,adultAges,e,0,R0,Cij,Ni,r,0*e,0*e,0*e,1,betaVaccine,0,1,1),spreadersAllocation,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
    optimalAllocationRes(ix)=computeFinalSize_generalized(optimalAllocation{ix},adultAges,e,0,R0,Cij,Ni,r,0*e,0*e,0*e,1,betaVaccine,0,1,1);
    if exitflag<=0
        optimalAllocation{ix}=spreadersAllocation;
        optimalAllocationRes(ix)=Inf;
    end

    if ix>1 
        [x,overallInfections_sample,exitflag,output]=fmincon(@(x) computeFinalSize_generalized(x,adultAges,e,0,R0,Cij,Ni,r,0*e,0*e,0*e,1,betaVaccine,0,1,1),optimalAllocation{ix-1},problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
        overallInfections_sample=computeFinalSize_generalized(x,adultAges,e,0,R0,Cij,Ni,r,0*e,0*e,0*e,1,betaVaccine,0,1,1);
        if exitflag>0 & overallInfections_sample<optimalAllocationRes(ix)
            optimalAllocation{ix}=x;
            optimalAllocationRes(ix)=overallInfections_sample;
        end
    end

    for kx=1:200
        
        y=randomizePoint(optimalAllocation{ix},Ni,Ni,e,vaccinesLeftToDistibute,0.05+mod(kx,3)/10);
        [x,overallInfections_sample,exitflag,output]=fmincon(@(x)computeFinalSize_generalized(x,adultAges,e,0,R0,Cij,Ni,r,0*e,0*e,0*e,1,betaVaccine,0,1,1),y,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
        overallInfections_sample=computeFinalSize_generalized(x,adultAges,e,0,R0,Cij,Ni,r,0*e,0*e,0*e,1,betaVaccine,0,1,1);
        if exitflag<=0
            overallInfections_sample=Inf;
        end

        if overallInfections_sample<optimalAllocationRes(ix)
            optimalAllocation{ix}=x;
            optimalAllocationRes(ix)=overallInfections_sample;
            %kx
        end
    end

    if vacOfLessSusceptible(ix)<optimalAllocationRes(ix)
        optimalAllocation{ix}=riskFocusedAllocation;
        optimalAllocationRes(ix)=vacOfLessSusceptible(ix);
    end
end

% Sweep back
for ix=(M-1):-1:1

 [x,overallInfections_sample,exitflag,output]=fmincon(@(x)computeFinalSize_generalized(optimalAllocation{ix+1},adultAges,e,0,R0,Cij,Ni,r,0*e,0*e,0*e,1,betaVaccine,0,1,1),y,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);

 if overallInfections_sample<optimalAllocationRes(ix)
        optimalAllocation{ix}=x;
        optimalAllocationRes(ix)=overallInfections_sample;
        'Improved at sweep back'
 end

end

return

