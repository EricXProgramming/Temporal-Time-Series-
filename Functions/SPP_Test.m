%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPP_Test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by:
% Eric Lai
% M.S. Student, Statistics
% University of California - Irvine
% Department of Statistics
% ellai@uci.edu

% Please send bug reports, comments, or questions to Eric Lai.
% This code comes with no guarantee or warranty of any kind.
% Last modified 8-4-2015. 

%% Notes
% This function will fit our data with the Stationary Poisson Process Model
% as well as analyze its fit by utilizing a transformed time series
% analysis supplemented by the Kolmogoro-Smirnov (K-S) test. 

% Inputs: 
% t -> A vector containing the timestamps for each event occurance in 
%      datenum format.

% Outputs: 
% mu -> The estimated parameter \mu for the constant rate at which events 
%       are expected to occur.
% AIC -> The Akaike Information Criterion value. 
% TransformedTimes -> A vector containing the transformed times for the
%                     time series. 
% U -> A vector containing the realizations (U_k) for the time series. 
% KSTest -> A vector containing the p-value, test-statistic, and critical 
%           value for the results of the K-S test respectively.     
 
function [mu,AIC,TransformedTimes,U,KSTest]=SPP_Test(t)
N = length(t);
T = max(t);

%% Compute the maximized likelihood estimate for the intensity function of a
%% Stationary Poisson Process model. 
mu=N/T; 

%% Compute the Akike Information Criterion value. 
L=N*log(mu)-mu*T; 
AIC=2*(1-L); 

%% Compute the Transformed Time \tau_k
TransformedTimes= t.*mu; % Compute the Transformed Time   
 
%% K-S Test for the Realizations U_k = 1-e^{\tau_{k}-\tau_{k-1}}
U=1-exp(-(TransformedTimes(2:end)-TransformedTimes(1:end-1)));
TestCDF=makedist('Uniform','Lower',0,'Upper',1); 
[~,P_value,Test_Statistic,Critical_Value] = kstest(U,'CDF',TestCDF);
KSTest=[P_value,Test_Statistic,Critical_Value];   
end