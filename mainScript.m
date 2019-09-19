%% Leaky integrate and fire model of neural spike generation

clear all

%% part 1: questions 1 to 3

gm = 0.5;       % muS/mm2
cm = 10;        % nF/mm2
r = 0.04;       % mm

A = 4*pi*r^2;   % surface area of cell, mm2

C = A*cm;       % membrane capacitance, nF
G = A*gm;       % membrane conductance, muS 
R = (1/G);      % membrane resistance, megaohm

%% part 3: questions 4 to 7

%% LIF model parameters are stored in a structure
p.C = 0.2;          % capacitance in nF
p.R = 100;          % resitance in megaohm
p.dt = 0.01;        % integration time step in milliseconds
p.dur = 0.5;        % simulation duration in sec
p.Vthresh = -60;    % threshold in mV
p.EL = -70;         % leakage reversal potential in mV
p.Vreset = -70;     % reset voltage in mV
p.V0 = -70;         % initial condition in mV

Ij = 0.15;          % injected current in nA 

[Vm, spikes] = myLIF(p,Ij);     % run the LIF model

t_vector = 0:p.dt:(length(Vm)-1)*p.dt;  % vector of time values

figure(1)           % plot LIF model results
plot(t_vector,Vm);
xlabel('time (s)');
ylabel('V_m (mV)');

f_rate = sum(spikes)/p.dur;     % firing rate of neuron

ISI = diff(find(spikes)).*p.dt; % inter-spike intervals
frateISI = (1./ISI)*1000;       % firing rate from ISIs, converted to Hz

%% questions 8 to 10

Ij = 0.15;          % injected current in nA
p.tref = 3;         % refractory period in msec

[Vm, spikes] = myLIFref(p,Ij);  % run LIF model with refractory period

t_vector = 0:p.dt:(length(Vm)-1)*p.dt;  % vector of time values

figure(2)           % plot new LIF model results
subplot(2,1,1)
plot(t_vector,Vm);
xlabel('time (s)');
ylabel('V_m (mV)');
prettyfigure;
subplot(2,1,2)
plot(t_vector,Vm);
xlim([0 100])
ylim([-75 -65])
xlabel('time (s)');
ylabel('V_m (mV)');
title('zoomed version');
prettyfigure;

f_rateref = sum(spikes)/p.dur;  % firing rate of neuron

ISI = diff(find(spikes)).*p.dt; % inter-spike intervals
frateISIref = (1./ISI)*1000;    % firing rate from ISIs, converted to Hz

% firing rate is slower because of the refractory period

%% questions 11 to 13

p.dur = 1;          % simulation duration in seconds
p.tref = 3;         % refractory period in msec
I = 0:0.01:0.5;     % vector of currents
h = length(I);      % number of current values
simul_rates = zeros(1,h);   % vector to store firing rates

    for k = 1:h         % run model for each current and store results
        [~,spikes] = myLIFref(p,I(k));
        simul_rates(k) = sum(spikes)/p.dur;
    end

rate = fiAnalytic(I,p);     % analytic solution

figure(3)                   % plot simulation results & analytic solution
plot(I*1000, rate, 'b')
hold on
plot(I*1000, simul_rates, 'rs')
hold off
xlabel('Injected current (pA)')
ylabel('Firing Rate (Hz)')
prettyfigure;

% simulation results match the analytic solution. The threshold current for
% a spike is 100 pA or 0.1 nA (the units used here)

%% questions 14 and 15

p.dur = 1;          % simulation duration in seconds
p.tref = 3;         % refractory period in msec
I = 0:0.1:10;       % vector of currents
h = length(I);      % number of current values
simul_rates = zeros(1,h);   % vector to store firing rates

    for k = 1 : h       % run model for each current and store results
        [~,spikes] = myLIFref(p,I(k));
        simul_rates(k) = sum(spikes)/p.dur;
    end

rate = fiAnalytic(I,p);     % analytic solution

figure(4)                   % plot simulation results & analytic solution
plot(I*1000, rate, 'b')
hold on
plot(I*1000,simul_rates, 'rs')
plot(I*1000, (1000/p.tref)*ones(h,1), '--k')
hold off
xlabel('Injected current (pA)')
ylabel('Firing Rate (Hz)')
prettyfigure;

% maximum firing rate possible is 1/tref, so 333.33 Hz

%% questions 16 to 18

I = 0.2;                % current
p.dur = 10;             % simulation duration in seconds
p.tref = 3;             % refractory period in msec
snoise = 0:0.05:0.4;    % standard deviations of noisy input current  
h = length(snoise);     % number of sample noise levels
ISI = cell(1,h);        % cell array to store ISIs for different noise levels

    for k = 1 : h       % run model for each noise level 
        [~,spikes] = myLIFnoise(p,I, snoise(k));
        spike_times = find(spikes == 1);
        spike_times = spike_times.*p.dt;
        ISI{k} = diff(spike_times);
    end
      
l2 = ceil(max(cellfun(@max,ISI)));      % compute bins for histograms
bins = 0:0.5:l2;
% compute counts of ISIs for histograms
counts = cellfun(@(x) histc(x,bins), ISI, 'UniformOutput', false);

figure(5)                   % plot histograms in nine-panel figure window
for k = 1 : length(counts)
    subplot(3,3,k)
    bar(bins, counts{k}./sum(counts{k}), 'histc');
    xlim([0 l2+5])
    if k > 3
        ylim([0 0.5]);
    end
    prettyfigure;
    xlabel('ISI (msec)');
    ylabel('probability');
    title(['\sigma_{noise} = ', num2str(1000*snoise(k)), ' (pA)']);
end

% compute mean and standard deviation of ISI distributions
mu = cellfun(@mean, ISI);
stdev = cellfun(@std, ISI);

figure(6)    % plot mean and standard deviation as a function of noise level
subplot(1,2,1)
plot(1000*snoise, mu);
prettyfigure;
title(' mean ISI for different \sigma_{noise}')
xlabel('\sigma_{noise} (pA)');
ylabel('mean ISI (msec)')
ylim([0 20]);
subplot(1,2,2)
plot(1000*snoise, stdev);
prettyfigure;
title(' std. dev. ISI for different \sigma_{noise}')
xlabel('\sigma_{noise} (pA)');
ylabel(' \sigma_{ISI} (msec)');

%% questions 19 and 20

p.dur = 1;          % duration of similation in seconds
p.tref = 3;         % refractory period in msec
I = 0:0.01:0.5;     % vector of input currents
h = length(I);      % number of input currents
noise_rates = zeros(1,h);   % vector to store firing rates

    for k = 1 : h   % run model for multiple input currents
       [~,spikes] = myLIFnoise(p,I(k), 0.4);
       noise_rates(k) = sum(spikes)/p.dur;
    end

figure(8)       % plot new model results
plot(1000*I, noise_rates, 's--k');
prettyfigure;
xlabel('\sigma_{noise} (pA)')
ylabel('Firing Rate (Hz)')
title('f-I curve for noise injection')

% prediction shifts threshold toward top left and makes it less abrupt and 
% less clear where the threshold is