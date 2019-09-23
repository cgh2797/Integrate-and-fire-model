

C = 0.2;          % capacitance in nF
R = 100;          % resitance in megaohm
dt = 0.01;        % integration time step in milliseconds
dur = 0.3;        % simulation duration in sec
Vthresh = -60;    % threshold in mV
EL = -70;         % leakage reversal potential in mV
Vreset = -70;     % reset voltage in mV
V0 = -70;         % initial condition in mV

tref = 3;         % refractory period in msec
Ij = 0.15;          % injected current in nA 


dur = dur*1000;      % simulation duration in msec
niter = floor(dur/dt)+1;    % number of iterations
V = EL;                     % initial condition
spikes = zeros(1,niter);    % vector for storing binary spike train
Vm = zeros(1,niter);        % vector for storing Vm
Vm(1) = V0;          % assign initial condition to first element of Vm
t_vector = 0:dt:(length(Vm)-1)*dt;  % vector of time values
tcounter = tref;            % counter for how long we have been refracted, 
                                % it is not refracted when we start
taum = R*C;                 % time constant in msec

    for idx =  2 : niter    % loop for number of iterations
        % check if we are refracted
        if tcounter < tref
            V = Vreset; 
            tcounter = tcounter + dt;   % update refractory counter
        else     % we integrate as before
            dVdt =(1/taum) .* ((EL - V) + R * Ij);
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

plot(t_vector,Vm);
xlabel('time (ms)');
ylabel('V_m (mV)');