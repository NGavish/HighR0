function SI_Fig1_R0star(collectData)
%% This function computes and produces the graph of SI Figure 1. R0* as function of VE and VC.
% This function should be run from main.m which defines some dependencies to the code

% Fix the value of sigma
sigma=2; N1=0.5; N2=0.5;

% Define grid
epsilonVec=0:0.025:1;epsilonVec(1)=[]; % Remove epsilon=0 since R_0^* blows up at epsilon->0
v=linspace(0,min(N1,N2),20);           % Recall vaccine coverage v<min(N1,N2)
e=v*0+1;
[epsilon,v]=meshgrid(epsilonVec,v);

%% Compute alpha(epsilon) 

% Function for alpha as defined in Theorem 2
myfun = @(x,sigma,epsilon) x-x^sigma-x^epsilon+x^(sigma*epsilon);

% For each value of epsilon, compute alpha(epsilon)
x0 = 0.5; % initial point
for ix=1:numel(epsilonVec)
    fun = @(x) myfun(x,sigma,epsilonVec(ix));
    alpha(ix)=fzero(fun,x0);
end

alpha=e'*alpha;

% Compute R0star
R0star=(N1+N2*sigma).*log(1./alpha)./(1-N1*alpha-N2*alpha.^sigma-v.*(alpha.^epsilon-alpha));

%% Presentation of the results
close all

% Tile layout
% subplot(2,2,1)

% Plot data
surf(epsilon,v,R0star);shading interp
hold on;
mesh(epsilon,v,R0star,'FaceColor','none','EdgeColor','k');
view(40,30)

% Window size
set(gcf,'Position',[961 396 474 401])

% Colors
colormap hsv

% Axes & title
xlabel('\epsilon');ylabel('v');zlabel('R_0^*    ','rotation',0);
set(gca,'xtick',[0.025 0.2:0.2:1]);axis tight
box on

printGraph('./graphs/R0star')

