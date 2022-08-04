%% Main function for producing all graphs of manuscript
% This function produces all graphs of the manuscript titled
% "Optimal vaccination at high reproductive numbers:sharp transitions and counterintuitive allocations"
% 
addpath('./AuxilaryFunctions/') 

% Output directories
[status,msg] = mkdir('./data');    % Computed data
[status,msg] = mkdir('./graphs');  % Graphs

collectData=false;

%% Figure 1. Example with three groups.
Fig1_ThreeGroupsExample(collectData)

%% Figure 2. Example using the USA demographic and social structure.
Fig2_USAExample(collectData)

%% Figure 3. Attack rate as a function of R_0 for different allocations.
Fig3_USAExample_attackRate(collectData)

%% Figure 4. Example with two groups.
Fig4_TwoGroupExample(collectData)

%% Figure S1. R0* as function of VE and VC.
SI_Fig1_R0star(collectData)

%% Figure S2. Rthreshold in various countries for different parameters.
SI_Fig2_USAExample(collectData)