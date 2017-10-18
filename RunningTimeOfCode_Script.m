%% Instructions
% This script file will compare the times of using two different methods to
% compute the transformed times (\tau_k) where
% \tau_k = \mu*t_k + \theta \sum_{t_i<t_k} \left( 1-e^{-\omega(t_k-t_i)} \right)
%        = \mu*t_k + \theta*(k-1)-\theta * e^{-\omega*(t_{k}-t_{k-1})}\left(
%          \sum_{t_{i}<t_{k-1}}e^{-\omega*(t_{k-1}-t_{i})}+1 \right) 
%% and Part1 =. 

%% Input values.
Threshold=0.05; % Input the desired threshold for a Tweet to be part of a topic. 
W = W; % Input the Tweet by Topic matrix that will be normalized. 
Times = datenum(times); % Input the time series in Day-Month-Year 
                        % Hour:Minute:Second (datetime) format.
Topic = 1; % Input the topic number.
Parameters=[1; 2; 3]; % Input the test parameters. 
                        
%% Normalize the rows of the Tweet by topic matrix.
NormalizedW = RowNormalizer(W);

%% Initialize the Time Stamps (note that the time stamps for Tweets may not be ordered).
Times = Times - Times(1); 

%% Compute a vector t containing times that correspond to a non-zero value 
%% in the chosen column of W above the Threshold. 
t=Times(NormalizedW(:,Topic)>Threshold);
t=sort(t);
N = length(t);
T = max(t);

tic; 

%% A faster method to calculate the transformed times but requires much more memory. 
% First consider the case that 1 < k \leq N.
tk = t(:,ones(N,1));
ti = t';
ti = ti(ones(N,1),:);
tmat = tk - ti;
tmat(1,:) = [];
Onemat = ones(size(tmat));
Summat = Parameters(2)*(Onemat-exp(-1*Parameters(3)*tmat));
Summat = tril(Summat,0);
TransformedTimes = Parameters(1)*t(2:end) + sum(Summat,2);
% Now add in the first case where k = 1.
TransformedTimes = [Parameters(1)*t(1); TransformedTimes(1:end)]; 

TransformedTimeM1_RunningTime=toc;
TransformedTimesM1 = TransformedTimes;

tic;

%% A slower method to calculate the transformed times but requires less memory. 
Sum = zeros(N,1); 
for k=2:N
    Sum(k) = exp(-Parameters(3)*(t(k)-t(k-1)))*(Sum(k-1)+1);
end 
Increments = (0:(N-1))';
TransformedTimes = Parameters(1)*t + Parameters(2)*(Increments-Sum);      

TransformedTimeM2_RunningTime=toc;
TransformedTimesM2 = TransformedTimes;

tic;

%% A faster method to calculate Part1 but requires much more memory.
% Case 1: k = 1
Part1 = log(Parameters(1)); 
% Case 2: 1 < k \leq N
tk = t(:,ones(N,1));
ti = t';
ti = ti(ones(N,1),:);
tmat = tk - ti;
tmat(1,:) = [];
Sum1mat = tril(Parameters(2)*Parameters(3)*exp(-1*Parameters(3)*tmat),0);
Sum1 = log(Parameters(1)+sum(Sum1mat,2));
Part1 = Part1+sum(Sum1);

Part1M1_RunningTime=toc;
Part1M1 = Part1;

tic;

%% A slower method to calculate Part1 but requires less memory.
Sum1 = zeros(N,1); 
for k=2:N
    Sum1(k) = exp(-Parameters(3)*(t(k)-t(k-1)))*(Sum1(k-1)+1);
end 
Part1=sum(log(Parameters(1)+Parameters(2)*Parameters(3)*Sum1));

Part1M2_RunningTime = toc;
Part1M2 = Part1;
