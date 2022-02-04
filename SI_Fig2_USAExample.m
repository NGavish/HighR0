function SI_Fig2_USAExample(collectData)
%% This function computes and produces the graph of SI Figure 2. Rthreshold in various countries for different parameters.
% This function should be run from main.m which defines some dependencies to the code

%% Computation of the optimal allocations
if collectData
    countryList={"BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"};

    % Construct list of combinations of parameters
    parameters=combvec(1:numel(countryList),55,[0.1 0.2 0.3]);
    parameters=[parameters combvec(1:numel(countryList),[40 50 55 60 70],0.2)];

    % For each combination, compute optimal allocations
    L=numel(parameters)/3;
    Rvalues=linspace(3,15,200);
    parfor ix=1:L
        countryName=countryList{parameters(1,ix)};
        vaccineCoverage=parameters(2,ix);
        varepsilon=parameters(3,ix);

        countryData=load(join(['./countryData/',countryName,'_data.mat'],''));

        Cij=countryData.contactMatrix;
        N=countryData.N;
        Ni=N*countryData.agDist;

        [vacOfMoreSusceptible,vacOfLessSusceptible,optimalAllocationRes,optimalAllocation,attackRatevacOfMoreSusceptible,attackRatevacOfLessSusceptible]=ComputeFinalSizeasFunctionOfR(Rvalues,Cij,Ni,vaccineCoverage,varepsilon);
        saveFile(join(['./data/data_',countryName,'_coverage=',num2str(vaccineCoverage),'_VE=',num2str(100*(1-varepsilon)),'.mat'],''),vacOfLessSusceptible,optimalAllocationRes,optimalAllocation,Rvalues,Ni,N,Cij);
    end

end

%% Presentation of the results
close all;

% Color scheme
defineColors;

% Load data
vaccineCoverage=55;

[countryList,countryListNames,Rthreshold01]=calcRthreshold(vaccineCoverage,0.1); % VE=90%
[countryList,countryListNames,Rthreshold]=calcRthreshold(vaccineCoverage,0.2); % VE=80%
[countryList,countryListNames,Rthreshold03]=calcRthreshold(vaccineCoverage,0.3); % VE=70%

varepsilon=0.2;

[countryList,countryListNames,Rthreshold40]=calcRthreshold(40,varepsilon); % VC=40%
[countryList,countryListNames,Rthreshold50]=calcRthreshold(50,varepsilon); % VC=50%
[countryList,countryListNames,Rthreshold60]=calcRthreshold(60,varepsilon); % VC=60%
[countryList,countryListNames,Rthreshold70]=calcRthreshold(70,varepsilon); % VC=70%

% Tile layout
tiledlayout(2,1);

%% First panel - Rthreshold for different VE
nexttile

% Plot data
plot(1:numel(countryList),Rthreshold01,'^',color=green,LineWidth=1);hold on;
plot(1:numel(countryList),Rthreshold,'o',color=blue,LineWidth=1)
plot(1:numel(countryList),Rthreshold03,'rs',LineWidth=1)

% Window size
set(gcf,'Position',[520 386 666 411])

% Axes & title
set(gca,'xticklabel',countryListNames);
axis([0.5 9.5 1 10.5]);grid on
ylabel('R_{threshold}')
title('R_{threshold} in various countries for different vaccine efficacies (VE)')

% Legend
legend('90% vaccine efficacy (VE)','80% VE','70% VE','location','south','Orientation','horizontal')

% Text
text(0.02,0.9,'A',Units='normalized',FontSize=13)

%% Second panel - Rthreshold for different VC
nexttile

% Plot data
plot(1:numel(countryList),Rthreshold70,'^',color=green,LineWidth=1);hold on;
plot(1:numel(countryList),Rthreshold60,'d',color=yellow,LineWidth=1);
plot(1:numel(countryList),Rthreshold,'o',color=blue,LineWidth=1)
plot(1:numel(countryList),Rthreshold50,'.',LineWidth=1);
plot(1:numel(countryList),Rthreshold40,'rs',LineWidth=1);

% Axes & title
set(gca,'xticklabel',countryListNames);
axis([0.5 9.5 1 13]);grid on
ylabel('R_{threshold}')
title('R_{threshold} in various countries for different vaccine coverage (VC)')

% Legend
legend('70% vaccine coverage (VC)','60% VC','55% VC','50% VC','40% VC','location','south','Orientation','horizontal')

% Text
text(0.02,0.9,'B',Units='normalized',FontSize=13)

% Export graphics
printGraph('./graphs/RthresholdVariousCountries')

function [countryList,countryListNames,Rthreshold]=calcRthreshold(vaccineCoverage,betaVaccine)
countryList={"BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"};
countryListNames={"Belgium", "USA", "India", "Spain", "Zimbabwe", "Brazil", "China", "South Africa", "Poland"};

for ix=1:numel(countryList)
    countryName=countryList{ix}
    data=load(join(['./data/data_',countryName,'_coverage=',num2str(vaccineCoverage),'_VE=',num2str(100*(1-betaVaccine)),'.mat'],''));
    Rthreshold(ix)=ComputeRthreshold(data);
end
return

function Rthreshold=ComputeRthreshold(data)
E=100*(1-(data.N-data.vacOfLessSusceptible)./(data.N-data.optimalAllocationRes));
%E=100*(1-data.optimalAllocationRes./data.vacOfLessSusceptible);
R0=data.Rvalues;
idx=1+max(find(E>1));
if idx<numel(R0)
    Rthreshold=R0(idx);
else
    Rthreshold=nan;
end
return