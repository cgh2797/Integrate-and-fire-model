function [Vm, spikes] = myLIFref(params, Iinj)
% [Vm, spikes] = myLIFref(params, Iinj)
% implements leaky integrate and fire model with a refractory period, given 
% input model parameters stored in a structure, and input current Iinj

C = params.C;               % capacitance in nF
R = params.R;               % resitance in megaohm
dt = params.dt;             % integration time step in msec
dur = params.dur*1000;      % simulation duration in msec
Vthresh = params.Vthresh;   % threshold in mV
EL = params.EL;             % leakage reversal potential in mV
Vreset = params.Vreset;     % reset voltage in mV
tref = params.tref;         % refractory period in msec
niter = floor(dur/dt)+1;    % number of iterations
V = EL;                     % initial condition
spikes = zeros(1,niter);    % vector for storing binary spike train
Vm = zeros(1,niter);        % vector for storing Vm
Vm(1) = params.V0;          % assign initial condition to first element of Vm

tcounter = tref;            % counter for how long we have been refracted, 
                                % it is not refracted when we start
taum = R*C;                 % time constant in msec

    for idx =  2 : niter    % loop for number of iterations
        % check if we are refracted
        if tcounter < tref
            V = Vreset; 
            tcounter = tcounter + dt;   % update refractory counter
        else     % we integrate as before
            dVdt =(1/taum) .* ((EL - V) + R * Iinj);
            V = V + dt .* dVdt;
        end 
        % check if spiking
        if V > Vthresh
            spikes(idx) = 1;
            V = Vreset;
            tcounter = 0;   % reset refractory counter
        end
        Vm(idx) = V;
    end

end