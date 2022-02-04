function Fig1_ThreeGroupsExample(collectData)
%% This function computes and produces the graph of Figure 1. Example with three groups.
% This function should be run from main.m which defines some dependencies to the code

%% Computation of the optimal allocations  
if collectData
    % Problem parameters
    Ni=[0.25 0.5 0.25]';    % Group sizes
    susceptiblity=[1 2 4];  % Susceptiblity profile (\sigma_i)
    varepsilon=0.2;         % Vaccine efficacy = 1-\varepsilon
    vaccineCoverage=40;     % Vaccine coverage

    % Define contact matrix
    Cij=diag(susceptiblity)*ones(size(Ni))*Ni'*diag(susceptiblity);
    [V,d]=eig(Cij);Cij=Cij/max(abs(diag(d)))

    % Vector of R0 values
    d=linspace(-0.1,0.1,51);
    Rvalues=unique([linspace(2.5,10,100) 4.72+d 5.82+d]);

    % Compute optimal vaccine allocations for vector of R0 values
    [vacOfMoreSusceptible,vacOfLessSusceptible,optimalAllocationRes,optimalAllocation,attackRatevacOfMoreSusceptible,attackRatevacOfLessSusceptible]=ComputeFinalSizeasFunctionOfR(Rvalues,Cij,Ni,vaccineCoverage,varepsilon);
    save './data/data3Groups,coverage=40,VE=80.mat'
end

%% Presentation of the results
close all;

% Color scheme
defineColors;   
newcolors=flipud(parula(4));
newcolors(3,:)=newcolors(1,:);
newcolors(1,:)=newcolors(2,:);
newcolors(2,:)=newcolors(3,:);
newcolors(3,:)=[];
colororder(newcolors)

% Load data and organize it in a matrix 
load './data/data3Groups,coverage=40,VE=80.mat'
Mat=flipud(cell2mat(optimalAllocation));

% Plot data
hndl=area(Rvalues,100*(Mat.*Ni)'/sum(Ni)) ;  

% Window size
set(gcf,'Position',[520 551 560 246])

% Colors
hndl(1).FaceAlpha=0.5;    hndl(2).FaceAlpha=0.5;

% Axes & title
axis([2.5 8.5 0 40])
xlabel('R_0');    
ylabel('allocation of vaccines');ytickformat('percentage');
title('Infections miniziming vaccine allocation as function of R_0')
set(gca,'ytick',0:10:40)
box on;

% Added text
text('Units','normalized','HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'String',{'High number of','daily interactions','(Group 3)'},...
    'Position',[0.19 0.3 0]);
text('Units','normalized','FontWeight','bold','FontSize',12,...
    'String','Group 2',...
    'Position',[0.40 0.5 0]);
text('Units','normalized','HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'String',{'Low number of','daily interactions','(Group 1)'},...
    'Position',[0.78 0.73 0],...
    'Color',[1 1 1]);

% Export graphics
printGraph(['./graphs/Three_group_example_coverage',num2str(vaccineCoverage)])
