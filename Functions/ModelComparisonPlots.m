%%%%%%%%%%%%%%%%%%%%%%%%% ModelComparisonPlots.m %%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by:
% Eric Lai
% M.S. Student, Statistics
% University of California - Irvine
% Department of Statistics
% ellai@uci.edu

% Please send bug reports, comments, or questions to Eric Lai.
% This code comes with no guarantee or warranty of any kind.
% Last modified 8-7-2015. 

% Inputs:
% EHP_Parameters -> A 1 x 3 vector containing the estimated parameters \mu, 
%                   \theta, and \omega respectively for an Exponential 
%                   Hawkes Process model.
% SPP_Parameters -> A 1 x 1 vector containing the estimated parameter \mu 
%                   for a Stationary Poisson process model.
% EHP_AIC -> The Akaike Information Criterion value for fitting an 
%            Exponential Hawkes Process model. 
% SPP_AIC -> The Akaike Information Criterion value for fitting an
%            Stationary Poisson Process model. 
% H -> A vector containing the timestamps for each event occurance in 
%      datenum format.
% Topic -> The index of the topic of consideration.

% Outputs: 
% A plot with two subplots containing estimated time series plots for the 
% Stationary Poisson Process and Exponential Hawkes Process model.

function ModelComparisonPlots(SPP_Parameter,EHP_Parameters,SPP_AIC,EHP_AIC,H,Topic)
%% Approximate times series plot for the Stationary Poisson Process model.
subplot(1,2,1);
plot([0 365], [SPP_Parameter SPP_Parameter],'b-','LineWidth',2);
xlabel('Time (in Days)','Fontsize',18); 
ylabel('Tweets/Day','Fontsize',18);
t=title(sprintf('Stationary Poisson Process with \n $\\hat{\\mu}$ = %.3e and AIC = %.3e',SPP_Parameter,SPP_AIC),'Fontsize',18);
set(t,'Interpreter','Latex','Fontsize',18);
axis([1,365,SPP_Parameter-1,SPP_Parameter+1]);

hold on;

%% Approximate times series plot for the Exponential Hawkes Process model. 
subplot(1,2,2);
Times = linspace(0,365,length(H));
N = length(Times); 
Sum = zeros(N,1);
for k=2:N
    Sum(k)=sum(exp(-EHP_Parameters(3)*(Times(k)-H(H<Times(k))))); 
end 
EHP_Points = EHP_Parameters(1) + EHP_Parameters(2)*EHP_Parameters(3)*Sum; 
plot(Times,EHP_Points,'b','LineWidth',2);
xlabel('Time (in Days)','Fontsize',18); 
ylabel('Tweets/Day','Fontsize',18);
t = title(sprintf('Exponential Hawkes Process with $\\hat{\\mu}$ = %.3e, \n $\\hat{\\theta}$ = %.3e, $\\hat{\\omega}$ = %.3e, and AIC = %.3e',EHP_Parameters(1),EHP_Parameters(2),EHP_Parameters(3),EHP_AIC),'Fontsize',18);
set(t,'Interpreter','Latex','Fontsize',18);
axis([1,365,min(EHP_Points)-1,max(EHP_Points)+1]);

%% Add a main title to the entire figure. 
t = suptitle(sprintf('Stationary Poisson Process vs. Exponential Hawkes Process for Topic %d',Topic));
set(t,'Interpreter','Latex','Fontsize',18);

saveas(gcf, strcat('EHP_Results/EHP_CP/CP_Topic_',num2str(Topic)), 'tiff'); 
saveas(gcf, strcat('EHP_Results/EHP_CP/CP_Topic_',num2str(Topic)), 'fig');
saveas(gcf, strcat('EHP_Results/EHP_CP/CP_Topic_',num2str(Topic)), 'epsc'); 
saveas(gcf, strcat('EHP_Results/EHP_CP/CP_Topic_',num2str(Topic)), 'fig'); 
end 