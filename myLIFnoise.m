function [Vm, spikes] = myLIFnoise(params, Iinj, noiseAmp)
% [Vm, spikes] = myLIFnoise(params, Iinj)
% implements leaky integrate and fire model with a refractory period and
% noisy input, given input model parameters stored in a structure, input 
% current Iinj, and input noise amplitude

C = params.C;               % capacitance in nF
R = params.R;               % resitance in megaohm
dt = params.dt;             % integration time step in msec
dur = params.dur*1000;      % simulation duration in msec
Vthresh = params.Vthresh;   % threshold in mV
EL = params.EL;             % leakage reversal potential in mV
Vreset= params.Vreset;      % reset voltage in mV
tref = params.tref;         % refractory period in msec
niter = floor(dur/dt)+1;    % number of iterations
V = EL;                     % initial condition
spikes = zeros(1,niter);    % vector for storing binary spike train
Vm = zeros(1,niter);        % vector for storing Vm
Vm(1)=params.V0;            % assign initial condition to first element of Vm

tcounter = tref;            % counter for how long we have been refracted, 
                                % it is not refracted when we start
taum = R*C;                 % time constant in msec            

% generate noisy input

    if noiseAmp == 0        % no noise
        input = Iinj*ones(1,niter);
    else                    % generate random noise
        noise = noiseAmp .* randn(1, niter);
        input = Iinj + noise;
    end

    for idx =  2 : niter    % loop for number of iterations
        % check if we are refracted
        if tcounter < tref
            V = Vreset;
            tcounter = tcounter + dt;
        else     % we integrate as before
            dVdt =(1/taum).* ( (EL - V)  + R * input(idx-1));
            V = V + dt .* dVdt;
        end
        % check if spiking
        if V > Vthresh
            spikes(idx) = 1;
            V = Vreset;
            tcounter = 0;   % reset refractory counter
        end 
        Vm(idx)=V;
    end

end