function [Vm, spikes] = myLIF(params, Iinj)
% [Vm, spikes] = myLIF(params, Iinj)
% implements leaky integrate and fire model, given input model parameters
% stored in a structure, and input current Iinj

C = params.C;               % capacitance in nF
R = params.R;               % resitance in megaohm
dt = params.dt;             % integration time step in msec
dur = params.dur*1000;      % simulation duration in msec
Vthresh = params.Vthresh;   % threshold in mV
EL = params.EL;             % leakage reversal potential in mV
Vreset = params.Vreset;     % reset voltage in mV
niter = floor(dur/dt)+1;    % number of iterations
V = EL;                     % initial condition
spikes = zeros(1,niter);    % vector for storing binary spike train
Vm = zeros(1,niter);        % vector for storing Vm 
Vm(1) = params.V0;          % assign initial condition to first element of Vm

taum = R*C;                 % time constant in msec
        
    for idx =  2 : niter    % loop for number of iterations
        % difference equation
        dVdt =(1/taum) .* ((EL - V) + R * Iinj);
        % update rule
        V = V + dt .* dVdt;
        % check if spiking
        if V > Vthresh
            spikes(idx) = 1;
            V = Vreset;
        end
        % update Vm
        Vm(idx) = V;
    end

end