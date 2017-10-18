%%%%%%%%%%%%%%%%%%%%%% EHP_ModelDiagnosticPlots.m %%%%%%%%%%%%%%%%%%%%%%%%%

% Written by:
% Eric Lai
% M.S. Student, Statistics
% University of California - Irvine
% Department of Statistics
% ellai@uci.edu

% Please send bug reports, comments, or questions to Eric Lai.
% This code comes with no guarantee or warranty of any kind.
% Last modified 8-6-2015. 

% Inputs:
% Parameters -> A vector containing the parameters \mu, \theta, and \omega
%               respectively for an Exponential Hawkes Process model.
% Topic -> The index of the topic of consideration.
% AIC -> The corresponding Akaike Information Criterion value. 
% TransformedTimes -> A vector containing the transformed times \tau_k 
%                     for the time series. 
% U -> A vector containing the realizations U_k for the time series. 

% Outputs: 
% A plot with four model diagnostic subplots.  

function EHP_ModelDiagnosticPlots(Parameters,Topic,AIC,TransformedTimes,U)
%% Model Diagnostic Plot 1
subplot(2,2,1);
Indices = 1:1:length(TransformedTimes);
plot(TransformedTimes,Indices,'k','LineWidth',2);
hold on;
subplot(2,2,1);
plot(Indices,Indices,'r','LineWidth',2);
x = xlabel('Transformed Times $\left( \tau_i \right)$','FontSize',18);
set(x,'Interpreter','Latex');
ylabel('Cumulative Number of Tweets','FontSize',18); 
axis([0,max(TransformedTimes),0,max(Indices)]);

%% Model Diagnostic Plot 2
N = length(U);
UCDF = [1:N]/N;
subplot(2,2,2);
plot(sort(U),UCDF,'k','LineWidth',2);
hold on;
subplot(2,2,2);
plot(UCDF,UCDF,'r','LineWidth',2);
x = xlabel('$U_{(k)}$','FontSize',18);
set(x,'Interpreter','Latex');
ylabel('Cumulative Distrubution','FontSize',18); 
axis([min(U),max(U),0,max(UCDF)]);

%% Model Diagnostic Plot 3
subplot(2,2,3);
hist(U)
x=xlabel('$U_{(k)}$','Fontsize',18); 
set(x,'Interpreter','Latex');
ylabel('Frequency','Fontsize',18) 

%% Model Diagnostic Plot 4
subplot(2,2,4);
hist3([U(1:N-1) U(2:N)],[50 50],'FaceAlpha',.5);
x=xlabel('$U_{(k)}$','Fontsize',18); 
set(x,'Interpreter','Latex');
y=ylabel('$U_{(k+1)}$','Fontsize',18); 
set(y,'Interpreter','Latex');
zlabel('Frequency','Fontsize',18); 

%% Add a main title to the entire figure. 
t = suptitle(sprintf('Topic %d Fitted with an Exponential Hawkes Process \n with $\\hat{\\mu}$ = %.3e, $\\hat{\\theta}$ = %.3e, $\\hat{\\omega}$ = %.3e, \n and AIC = %.3e',Topic,Parameters(1),Parameters(2),Parameters(3),AIC));
set(t,'Interpreter','Latex','Fontsize',18);

saveas(gcf, strcat('EHP_Results/EHP_MDP/EHP_MDP_Topic_',num2str(Topic)), 'tiff'); 
saveas(gcf, strcat('EHP_Results/EHP_MDP/EHP_MDP_Topic_',num2str(Topic)), 'epsc'); 
saveas(gcf, strcat('EHP_Results/EHP_MDP/EHP_MDP_Topic_',num2str(Topic)), 'fig'); 