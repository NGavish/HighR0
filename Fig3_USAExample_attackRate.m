function Fig3_USAExample_attackrate(collectData)
%% This function computes and produces the graph of Figure 3. Attack rate as a function of R_0 for different allocations.
% This function should be run from main.m which defines some dependencies to the code

%% Computation of the optimal allocations
if collectData
        varepsilon=0.2;                % Vaccine efficacy = 1-\varepsilon
        vaccineCoverage=55;            % Vaccine coverage

    if ~isfile('./data/dataUS9by9,coverage=55,VE=80.mat')
        % Load USA data
        countryData=load('./CountryData/USA_data.mat');

        % Problem parameters
        Cij=countryData.contactMatrix; % Contact matrix
        N=countryData.N;               % Total population
        Ni=N*countryData.agDist;       % Group sizes

        % Vector of R0 values
        Rvalues=unique([linspace(1,10,120) 4.9+linspace(-0.1,0.1,51) 5.53+linspace(-0.1,0.1,51) 5.9+linspace(-0.1,0.1,51) 7.46+linspace(-0.1,0.1,51)]);

        % Compute optimal vaccine allocations for vector of R0 values
        [vacOfMoreSusceptible,vacOfLessSusceptible,optimalAllocationRes,optimalAllocation,attackRatevacOfMoreSusceptible,attackRatevacOfLessSusceptible]=ComputeFinalSizeasFunctionOfR(Rvalues,Cij,Ni,vaccineCoverage,varepsilon);
        saveFile(['./data/dataUS9by9_coverage=',num2str(vaccineCoverage),',VE=',num2str(100*(1-varepsilon)),'.mat'],vacOfLessSusceptible,optimalAllocationRes,optimalAllocation,Rvalues,Ni,N,Cij);
    end

    % Load & analyze data
    data=load('./data/dataUS9by9,coverage=55,VE=80.mat');

    % Extract demographic parameters
    Ni=data.Ni;
    N=numel(Ni);
    populationSize=sum(Ni);
    adultAges=1:N;
    e=ones(size(Ni));

    % Extract `spreaders' allocation
    % This is the optimal allocation at smallest R0 in which herd immunity is not reached
    idx=22; % This is the index of the R0\approx 2.6
    optimalAllocationAtLowR=data.optimalAllocation{idx};

    % Extract final sizes resulting from optimal and asymptotic
    % allocations
    R0values=data.Rvalues(idx:end);
    finalSize_optimal=data.optimalAllocationRes(idx:end);
    finalSize_asymptotic=data.vacOfLessSusceptible(idx:end);

    % Compute final sizes resulting from spreaders and uniform allocations
    for ix=1:numel(R0values)
        R0=R0values(ix);
        [finalSize_spreaders(ix),dummy,dummy,dummy]=computeFinalSize_generalized(optimalAllocationAtLowR,adultAges,e,0,R0,data.Cij,Ni,0*e,0*e,0*e,0*e,1,varepsilon,0,1,1);
        [finalSize_uniform(ix),dummy,dummy,dummy]=computeFinalSize_generalized(vaccineCoverage*e/100,adultAges,e,0,R0,data.Cij,Ni,0*e,0*e,0*e,0*e,1,varepsilon,0,1,1);
    end
    N=populationSize;
    save('./data/USAExample_attackRate','R0values','finalSize_spreaders','finalSize_uniform','finalSize_asymptotic','finalSize_optimal','N')
end

%% Presentation of the results
close all;
% Color scheme
defineColors;

% Load data and organize it in a matrix
data=load('./data/USAExample_attackRate');

% Plot data
plot(data.R0values,100*data.finalSize_spreaders/data.N,color=blue,LineStyle='-.',LineWidth=1.5);hold on;
plot(data.R0values,100*data.finalSize_asymptotic/data.N,Color=red,LineWidth=1.5);
plot(data.R0values,100*data.finalSize_optimal/data.N,'k--',linewidth=1.5);hold on;
plot(data.R0values,100*data.finalSize_uniform/data.N,':',Color=green,linewidth=1.5)

% Window size
set(gcf,'Position',[520 551 560 246])

% Axes & title
xlabel('R_0');ylabel('attack rate')
ylim([0 100]);
ytickformat('percentage');
title('Attack rate - USA example')
box on;grid on;
set(gca,'xtick',[data.R0values(1) 3:10],'xticklabel',{'2.6','3','4','5','6','7','8','9','10'})
axis([2.5 10.1 0 100]);

% Legend
legend('Vaccination of spreaders','Asymptotic allocation','Optimal vaccine allocation','Uniform vaccination','location','northwest')

% Export graphics
printGraph(['./graphs/USExample_attackrate'])
