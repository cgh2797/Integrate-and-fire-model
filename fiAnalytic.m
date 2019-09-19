function frate = fiAnalytic(I, params)
% frate = fiAnalytic(I, params)
% implements analytic solution for f-I curve for leaky integrate and fire
% model of neural spike generation, given input vector of currents and
% input structure of model parameters

Vthresh = params.Vthresh;       % threshold in mV
Vreset = params.Vreset;         % reset voltage in mV
gL = 1/params.R;                % membrane conductance
C = params.C;                   % capacitance in nF
taum = C/gL;
frate = zeros(1, length(I));    % firing rate
rheobase = gL*(Vthresh-Vreset);
tref = params.tref;             % refractory period in msec
frate(I>rheobase) = 1000./(tref+taum .* log(I(I>rheobase)./(I(I>rheobase) - gL .* (Vthresh-Vreset))));

end