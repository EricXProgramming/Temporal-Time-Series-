%%%%%%%%%%%%%%%%%%%%%%% Time_Series_Script.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by:
% Eric Lai
% M.S. Student, Statistics
% University of California - Irvine
% Department of Statistics
% ellai@uci.edu

% Please send bug reports, comments, or questions to Eric Lai.
% This code comes with no guarantee or warranty of any kind.
% Last modified 8-6-2015. 

%% Notes 
% Run this script file in order to test whether or not our data 
% fits the Stationary Poisson Process and Exponential Hawkes Process well.
% The Akaike Information Criterion value, results from the 
% KS-Test test, parameter estimates, model diagnostic plots, and multiple
% comparison plots between each point process model will be generated. 
% Topics with only no Tweets or one Tweet are excluded from this analysis. 
% The transformed times \tau_k and and realizations U_k that have values of
% infinity will not be saved. 

% Abbreviations: 
% SPP -> Stationary Poisson Process
% EHP -> Exponential Hawkes Process
% CP -> Comparison Plots
% AIC -> Akaike Information Criterion
% KS-Test -> Kolmogorov–Smirnov Test

% Parameters for the SPP/EHP Model:
% 1. \mu -> stationary rate of occurrences for background events.
%2. \theta -> mean number events triggered by an arbitrary event.
% 3. \omega^{-1} -> expected amount of time it takes for an event to 
%                   trigger another event.

% Variables that are Saved:
% NumberOfTweetsInTopic -> A N x 1 vector containing the number of Tweets
%                          in each topic for a given threshold where entry 
%                          (i,1) corresponds to the number of Tweets in 
%                          topic i. 
% EventOccurances -> A N x 1 cell where the (i,1) entry corresponds to a 
%                    vector containing the times that correspond to a 
%                    non-zero value in column i (topic i) of W above the 
%                    Threshold.  
% SPP_Parameter -> A N x 1 vector containing the estimated parameter for a
%                  SPP model. Entry (i,1) gives the estimated parameter \mu
%                  for topic i.
% SPP_AIC -> A N x 1 vector where entry (1,i) gives the Akaike Information 
%            Criterion value for fitting topic i with an SPP model. 
% SPP_KSTest -> A N x 3 matrix containing the results of the K-S Test to
%               determine the fit of an SPP model. The entries (i,1), 
%               (i,2), and (i,3) are the P-value, test statistic, and 
%               critical value respectively for the K-S Test and topic i.
% SPP_TransformedTimes -> A N x 1 cell where entry (i,1) contains the 
%                         transformed times \tau_k for fitting an SPP model 
%                         with topic i. 
% SPP_U_Topic -> A N x 1 cell where entry (i,1) gives the realizations 
%                U_k for fitting an SPP model for topic i.
% SPP_TimeToTest -> The time in hours it takes to test if our data fits an
%                   SPP sufficiently well.
% EHP_Parameters -> A N x 3 matrix containing the estimated parameters an 
%                   EHP model. Entries (i,1), (i,2), and (i,3) corresponds 
%                   to the estimated parameters \mu, \theta, and \omega for 
%                   topic i, respectively.
% EHP_AIC -> A N x 1 vector where entry (i,1) gives the Akaike Information 
%            Criterion value for fitting topic i with an EHP model.
% EHP_KSTest -> A N x 3 vector containing the results of the K-S Test to
%               determine the fit of an EHP model. The entries (i,1), 
%               (i,2), and (i,3) are the P-value, test statistic, and 
%               critical value respectively for the K-S Test and topic i. 
% EHP_TransformedTimes -> A N x 1 cell where entry (i,1) contains the 
%                         transformed times \tau_k for fitting an EHP 
%                         model with topic i. 
% EHP_U_Topic -> A N x 1 cell where entry (i,1) gives the realizations 
%                U_k for fitting an EHP model for topic i.
% EHP_TimeToTest -> The time in hours it takes to test if our data fits an
%                   EHP sufficiently well.

%% Input values for the script file. 
InitialTopic=1; % Input the Topic number to start at.
EndingTopic=300; % Input the Topic number to end at. 
Threshold=0.05; % Input the desired threshold for a Tweet to be part of a topic. 
W = W; % Input the Tweet by Topic matrix that will be normalized along the rows. 
Times = sort(datenum(times)); % Input the time series as a column vector in 
                              % Day-Month-Year Hour:Minute:Second 
                              %(datetime) format and in increasing order.
                        
%% Create some folders and subfolders to save the generated .mat and .fig 
%% files. 
mkdir('Time Series Results')
mkdir('SPP_Results','SPP_U')
mkdir('SPP_Results','SPP_TransformedTimes')
mkdir('SPP_Results','SPP_MDP') 
mkdir('SPP_Results','SPP_CP')
mkdir('EHP_Results','EHP_U')
mkdir('EHP_Results','EHP_TransformedTimes')
mkdir('EHP_Results','EHP_MDP')
mkdir('EHP_Results','EHP_CP')
                        
%% Normalize the rows of the Tweet by topic matrix
NormalizedW = RowNormalizer(W);

%% Initialize the Time Stamps.
Times = Times - Times(1); 

%% Compute the number of Tweets in each topic and gather the indices that 
%% do not correspond to a Topic with no Tweets or one Tweet and save it. 
%% Afterwards, compute a cell t_Cell containing times that correspond to a 
%% non-zero value in the chosen column of W above the Threshold and save it.  
AllIndices = InitialTopic:EndingTopic;
NumberOfTweetsInTopic = zeros(EndingTopic,1);
t_Cell = cell(EndingTopic,1);
for i=AllIndices
    t=Times(NormalizedW(:,i)>Threshold);
    NumberOfTweetsInTopic(i) = length(t);
    t_Cell{i}=t;
end     
BadIndices = find(NumberOfTweetsInTopic==0 | NumberOfTweetsInTopic==1);
Indices = AllIndices(~ismember(AllIndices,BadIndices));
save('Time Series Results/NumberOfTweetsInTopic.mat','NumberOfTweetsInTopic');
save('Time Series Results/EventOccurances.mat','t_Cell')

tic;

%% Test whether or not a SPP model is a good fit. 
SPP_Parameter = zeros(EndingTopic,1);
SPP_AIC = zeros(EndingTopic,1);
SPP_KSTest = zeros(EndingTopic,3);
SPP_Cell_U = cell(EndingTopic,1);
SPP_Cell_TransformedTimes = cell(EndingTopic,1);
for i=Indices   
    [SPP_Parameter(i),SPP_AIC(i),SPP_TransformedTimes,SPP_U,SPP_KSTest(i,:)]=SPP_Test(t_Cell{i});
    if sum(isinf(SPP_TransformedTimes))~=0 || sum(isinf(SPP_U))~=0
        continue
    end 
    SPP_ModelDiagnosticPlots(SPP_Parameter(i),i,SPP_AIC(i),SPP_TransformedTimes,SPP_U); 
    SPP_Cell_TransformedTimes{i} = SPP_TransformedTimes;
    SPP_Cell_U{i} = SPP_U;
end

SPP_TimeToTest = toc;

%% Save the results from fitting a SPP model as a .mat file. 
save('SPP_Results/SPP_Parameter.mat','SPP_Parameter');
save('SPP_Results/SPP_AIC.mat','SPP_AIC');
save('SPP_Results/SPP_KSTest.mat','SPP_KSTest');
save('SPP_Results/SPP_TimeToTest.mat','SPP_TimeToTest');
save('SPP_Results/SPP_TransformedTimes/SPP_TransformedTimes','SPP_Cell_TransformedTimes');
save('SPP_Results/SPP_U/SPP_U','SPP_Cell_U'); 

tic;

%% Test whether or not an EHP model is a good fit. 
EHP_Parameters = zeros(EndingTopic,3);
EHP_AIC = zeros(EndingTopic,1);
EHP_KSTest=zeros(EndingTopic,3);
EHP_Cell_U = cell(EndingTopic,1);
EHP_Cell_TransformedTimes = cell(EndingTopic,1);
for i=Indices
    [EHP_Parameters(i,:),EHP_AIC(i),EHP_TransformedTimes,EHP_U,EHP_KSTest(i,:)]=EHP_Test(t_Cell{i});
    if sum(isinf(EHP_TransformedTimes))~=0 || sum(isinf(EHP_U))~=0
        continue
    end 
    EHP_ModelDiagnosticPlots(EHP_Parameters(i,:),i,EHP_AIC(i),EHP_TransformedTimes,EHP_U);
    EHP_Cell_TransformedTimes{i} = EHP_TransformedTimes;
    EHP_Cell_U{i} = EHP_U;
end
 
EHP_TimeToTest = toc;
 
%% Save the results from fitting an EHP model as a .mat file. 
save('EHP_Results/EHP_Parameters.mat','EHP_Parameters');
save('EHP_Results/EHP_AIC.mat','EHP_AIC');
save('EHP_Results/EHP_KSTest.mat','EHP_KSTest');
save('EHP_Results/EHP_TimeToTest.mat','EHP_TimeToTest');
save('EHP_Results/EHP_TransformedTimes/EHP_TransformedTimes','EHP_Cell_TransformedTimes');
save('EHP_Results/EHP_U/EHP_U','EHP_Cell_U'); 

%% Create a comparison plot of the P-values between a SPP model and an EHP model. 
figure()
set(gcf, 'Color', 'w');
plot(SPP_KSTest(:,1),'g','LineWidth',2);
hold on
plot(EHP_KSTest(:,1),':.b','LineWidth',2);
hold on 
plot([0 300], [0.05 0.05],'r','LineWidth',2)
xlabel('Topic Number', 'Fontsize',18);
ylabel('P-value','Fontsize',18); 
legend('Stationary Poisson Process', 'Exponential Hawkes Process','P-value = 0.05','Location','northeast','Fontsize',18)
saveas(gcf,'Time Series Results/Pvalue_MC','tiff')
saveas(gcf,'Time Series Results/Pvalue_MC','epsc') % epsc is an eps file with color 
saveas(gcf,'Time Series Results/Pvalue_MC','fig')

%% Create a comparison plot of the AIC values between a SPP model and an EHP model. 
figure()
plot(SPP_AIC,'g','LineWidth',2);
hold on
plot(EHP_AIC,':.b','LineWidth',2);
xlabel('Topic Number', 'Fontsize',18);
ylabel('Akaike Information Criterion Value','Fontsize',18); 
legend('Stationary Poisson Process', 'Exponential Hawkes Process','Location','northeast','Fontsize',18)
saveas(gcf,'Time Series Results/AICvalue_MC','tiff')
saveas(gcf,'Time Series Results/AICvalue_MC','epsc') % epsc is an eps file with color 
saveas(gcf,'Time Series Results/AICvalue_MC','fig')

%% Create comparison plots involving the different time series models for each topic.
for i=Indices
    ModelComparisonPlots(SPP_Parameter(i),EHP_Parameters(i,:),SPP_AIC(i),EHP_AIC(i),t_Cell{i},i);
end 
    
