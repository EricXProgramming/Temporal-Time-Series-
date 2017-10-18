%%%%%%%%%%%%%%%%%%%%%%%%%% EHP_NegativeLogLikelihood.m %%%%%%%%%%%%%%%%%%%%
% Written by:
% Eric Lai
% M.S. Student, Statistics
% University of California - Irvine
% Department of Statistics
% ellai@uci.edu

% Please send bug reports, comments, or questions to Eric Lai.
% This code comes with no guarantee or warranty of any kind.
% Last modified 8-4-2015.

% This function will compute the negative of the log-likelihood function 
% for an Exponential Hawkes process model given the following input values. 

% Input: 
% Parameters -> A vector containing the parameters \mu, \theta, and \omega 
%               for an exponential Hawkes process model, respectively. 
% t -> The time series data in datenum format. 

% Output:
% LogLikelihood -> The value of the log-likelihood function that was
%                  computed. 

%% Notes: 
% This function will code the following: 
% \mu -> stationary rate of occurrences for background events.
%\theta -> mean number events triggered by an arbitrary event.
% \omega^{-1} -> expected amount of time it takes for an event to trigger 
%                another event.
% Part1 -> \sum_{k=1}^{N} log\left( \sum _{i=1}^{t_i<t_k} \theta*\omega e^{-\omega(t_k-t_i)} \right)
%          = \sum_{k=1}^N log \left( e^{-\omega*(t_{k}-t_{k-1})} \left(
%            \sum_{t_{i}<t_{k-1}}e^{-\omega*(t_{k-1}-t_{i})}+1 \right)
%            \right)
% Part2 -> \mu*T \theta*\sum _{i=1}^{N} \left(1-e^{-\omega(T-t_i)}\right)

function [NegativeLogLikelihood]=EHP_NegativeLogLikelihood(Parameters,t)
T=max(t);
N=length(t);

%% Calculating Part1 using a faster method that requires less memory.
Sum1 = zeros(N,1); 
for k=2:N
    Sum1(k) = exp(-Parameters(3)*(t(k)-t(k-1)))*(Sum1(k-1)+1);
end 
Part1=sum(log(Parameters(1)+Parameters(2)*Parameters(3)*Sum1));

%% Calculating Part2
Tvec = T*ones(size(t));
Sum2 = ones(size(t))-exp(-Parameters(3)*(Tvec-t));
Part2=Parameters(1)*T+Parameters(2)*sum(Sum2);

%% Calculating the negative of the log-likelihood function 
NegativeLogLikelihood=-(Part1-Part2);
end 



