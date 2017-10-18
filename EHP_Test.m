%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EHP_Test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% This function will fit our data with the Exponential Hawkes 
% Process Model as well as analyze its fit by utilizing a transformed Times 
% series analysis supplemented by the Kolmogoro-Smirnov (K-S) test. 

% Inputs: 
% t -> A vector containing the timestamps for each event occurance in 
%      datenum format.

% Outputs: 
% Parameters -> A vector containing the parameters \mu, \theta, and \omega 
%               for an Exponential Hawkes Process model, respectively. 
% AIC -> The Akaike Information Criterion value for fitting an Exponential
%        Hawkes Process model. 
% TransformedTimes -> A vector containing the transformed times \tau_k 
%                      for the times series. 
% U -> A vector containing the realizations U_k for the times series. 
% KSTest -> A vector containing the p-value, test-statistic, and critical 
%           value for the results of the K-S test, respectively.     

function [Parameters,AIC,TransformedTimes,U,KSTest]=EHP_Test(t)
N = length(t);
T = max(t);

%% Create and modify the optimization problem. 
x0=[N/(2*T),0.5,1];

% Uncomment the following to lines to use fmincon as the optimization algorithm. 
fminconOptions = optimoptions('fmincon','Algorithm','interior-point','Diagnostics','off','UseParallel','never');
[Parameters,MinLogLikelihood] = fmincon(@(Parameters)EHP_NegativeLogLikelihood(Parameters,t),x0,[],[],[],[],[0,0,0],[100,100,100],[],fminconOptions); 

% Uncomment the following two lines to use fminunc as the optimization
% algorithm. Note that fminunc may give rise to parameters that are
% incredibly large, since this is an unconstrained optimization. 
% fminuncOptions = optimoptions('fminunc','Algorithm','quasi-newton','Diagnostics','off');
% [Parameters,MinLogLikelihood]=fminunc(@(Parameters)EHP_NegativeLogLikelihood(Parameters,t),x0,fminuncOptions);

MaxLogLikelihood=-MinLogLikelihood;
AIC=6-2*MaxLogLikelihood;

%% Compute the transformed times 
%% \tau_k = \mu \cdot t_k + \theta \sum_{t_i<t_k} \left( 1-
%%          e^{-\omega(t_k-t_i)} \right)
%%        = \mu \cdot t_k + \theta \cdot (k-1)-\theta \cdot e^{-\omega 
%%          \cdot (t_{k}-t_{k-1})}\left(
%%          \sum_{t_{i}<t_{k-1}}e^{-\omega \cdot (t_{k-1}-t_{i})}+1 \right) 
%% using a faster method that requires less memory.
Sum = zeros(N,1); 
for k=2:N
    Sum(k) = exp(-Parameters(3)*(t(k)-t(k-1)))*(Sum(k-1)+1);
end 
Increments = (0:(N-1))';
TransformedTimes = Parameters(1)*t + Parameters(2)*(Increments-Sum);

%% K-S Test for the Realizations U_k = 1-e^{\tau_{k}-\tau_{k-1}}
U=1-exp(-(TransformedTimes(2:end)-TransformedTimes(1:end-1)));
TestCDF=makedist('Uniform','Lower',0,'Upper',1); 
[~,P_value,Test_Statistic,Critical_Value] = kstest(U,'CDF',TestCDF);
KSTest=[P_value,Test_Statistic,Critical_Value];   
end 