%%%%%%%%%%%%%%%%%%%%%%%%%% SPP_Simulation.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eventTimes = SPP_Simulation(lambda,T)
% This function returns a column vector of event times for a simulated 
% stationary poisson process over observation window [0,T] with constant rate 
% lambda
% Inputs:
% lambda        = the rate of the poisson process
% T             = time over which to simulate
%
numEvents = poissrnd(lambda*T);
uniforms = rand(numEvents,1);
eventTimes = T * uniforms;
eventTimes = sort(eventTimes);
end 

    
    


