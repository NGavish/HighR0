function Fig2_USAExample(collectData)
%% This function computes and produces the graph of Figure 2. Example using the USA demographic and social structure.
% This function should be run from main.m which defines some dependencies to the code

%% Computation of the optimal allocations
if collectData
    % Load USA data
    countryData=load('./CountryData/USA_data.mat');

    % Problem parameters
    Cij=countryData.contactMatrix; % Contact matrix
    N=countryData.N;               % Total population
    Ni=N*countryData.agDist;       % Group sizes

    varepsilon=0.2;                % Vaccine efficacy = 1-\varepsilon
    vaccineCoverage=55;            % Vaccine coverage

    % Vector of R0 values
    Rvalues=unique([linspace(1,10,120) 4.9+linspace(-0.1,0.1,51) 5.53+linspace(-0.1,0.1,51) 5.9+linspace(-0.1,0.1,51) 7.46+linspace(-0.1,0.1,51)]);

    % Compute optimal vaccine allocations for vector of R0 values
    [vacOfMoreSusceptible,vacOfLessSusceptible,optimalAllocationRes,optimalAllocation,attackRatevacOfMoreSusceptible,attackRatevacOfLessSusceptible]=ComputeFinalSizeasFunctionOfR(Rvalues,Cij,Ni,vaccineCoverage,varepsilon);
    saveFile(['./data/dataUS9by9_coverage=',num2str(vaccineCoverage),',VE=',num2str(100*(1-varepsilon)),'.mat'],vacOfLessSusceptible,optimalAllocationRes,optimalAllocation,Rvalues,Ni,N,Cij);
end

%% Presentation of the results
close all;

% Color scheme
defineColors;

% Load data and organize it in a matrix
data=load('./data/dataUS9by9,coverage=55,VE=80.mat');
Mat=cell2mat(data.optimalAllocation);

% Plot data
idx=find(data.optimalAllocationRes>1e-4*data.N);
hndl=area(data.Rvalues(idx),100*(Mat(:,idx).*data.Ni)'/data.N);

% Window size
set(gcf,'Position',[520 551 560 246])

% Colors
for ix=1:9
    hndl(ix).FaceColor=lineColors(ix,:);
end

% Axes & title
axis tight
xlabel('R_0');ylabel('allocation of vaccines')
ytickformat('percentage')
title('Infections minimizing vaccine allocations - USA example')
set(gca,'xtick',[data.Rvalues(22) 3:10],'xticklabel',{'2.6','3','4','5','6','7','8','9','10'})
box on;

% Added text
textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};

% Add text at R0\approx 2.6
idx=22;
currY=0;
for ix=1:7
    currY=currY+100*Mat(ix,idx).*data.Ni(ix)/data.N;
    if Mat(ix,idx)>sum(Mat(:,idx))/60
        h(ix)=text(data.Rvalues(idx),currY,['  ',textAges{ix},' (',num2str(100*Mat(ix,idx),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold')
    end
end
hold on;plot([data.Rvalues(idx) data.Rvalues(idx)],[0 55],'w:','linewidth',2)

% Add text at R0\approx 8.6
idx=305;
currY=0;
for ix=1:9
    currY=currY+100*Mat(ix,idx).*data.Ni(ix)/data.N;
    if Mat(ix,idx)>sum(Mat(:,idx))/60
        h(ix)=text(data.Rvalues(idx),currY,['  ',textAges{ix},' (',num2str(100*Mat(ix,idx),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold')
    end
end
hold on;plot([data.Rvalues(idx) data.Rvalues(idx)],[0 55],'w:','linewidth',2)

% Export graphics
printGraph(['./graphs/USIntroExample'])
