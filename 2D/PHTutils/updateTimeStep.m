function [ nextTimeStep,time,previousTime] = updateTimeStep( timeStep,time,error)

disp('updating time step')
rho=0.9;
tol=1e-3;
%nextTimeStep=rho*(tol/error)^(0.5)*timeStep;

hfac_lte = rho*(tol/error)^(0.5);

%clip to interval [0.1, 10]
if hfac_lte<0.1
    hfac_lte = 0.1;
elseif hfac_lte>10
    hfac_lte = 10;
end

nextTimeStep = hfac_lte*timeStep;

%constrain time step
if nextTimeStep> 1e-6
    nextTimeStep=1e-6;
end

% if (time+timeStep+nextTimeStep >=max_time)
%     nextTimeStep=max_time - (time + timeStep);
% end

previousTime=time;
time=time+timeStep;

end

