function Fig4_TwoGroupsExample(collectData)
%% This function computes and produces the graph of Figure 4. Example with two groups.
% This function should be run from main.m which defines some dependencies to the code

%% Computation of the optimal allocations  
if collectData
    % Problem parameters
    Ni=[0.5 0.5]';    % Group sizes
    susceptiblity=[1 2];  % Susceptiblity profile (\sigma_i)
    varepsilon=0.2;         % Vaccine efficacy = 1-\varepsilon
    vaccineCoverage=40;     % Vaccine coverage

    % Define contact matrix
    Cij=diag(susceptiblity)*ones(size(Ni))*Ni';
    [V,d]=eig(Cij);Cij=Cij/max(abs(diag(d)));

    % Vector of R0 values
    Rvalues=unique([linspace(1,5,100) 1.8+linspace(-0.05,0.05,11) 3.3+linspace(-0.05,0.05,101)]);
    [vacOfMoreSusceptible,vacOfLessSusceptible,optimalAllocationRes,optimalAllocation,attackRatevacOfMoreSusceptible,attackRatevacOfLessSusceptible]=ComputeFinalSizeasFunctionOfR(Rvalues,Cij,Ni,vaccineCoverage,varepsilon);

    save './data/data2Groups,coverage=40,VE=80.mat'
end

%% Presentation of the results
close all;

% Color scheme
defineColors;   

% Load data and organize it in a matrix
data=load('./data/data2Groups,coverage=40,VE=80.mat');

% Plot data
plot(data.Rvalues,100*data.vacOfMoreSusceptible,color=blue,LineStyle='-.',LineWidth=1.5);hold on;
plot(data.Rvalues,100*data.vacOfLessSusceptible,Color=red,LineWidth=1.5);
plot(data.Rvalues,100*data.optimalAllocationRes,'k--',linewidth=1.5);hold on;

% Window size
set(gcf,'Position',[520 551 560 246])

% Axes & title
xlabel('R_0');ylabel('attack rate')
ylim([0 100]);
ytickformat('percentage');
title('attack rate - two groups example')
box on;grid on;
set(gca,'xtick',1:5)
axis([1 5 0 100]);

% Legend
legend('Vaccination of the more susceptible (spreaders)','Vaccination of the less susceptible (asymptotic allocation)','Optimal vaccine allocation','location','northwest')

% Export graphics
printGraph(['./graphs/TwoGroupExample'])
