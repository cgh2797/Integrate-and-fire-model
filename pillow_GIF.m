function [tsp,Vmem,Ispk,Istm] = gifrunmod(gifprs, Stim, dtsim);
% [tsp,Vmem,Ispk,Istm] = gifrunmod(gifprs, Stim, dtsim);
%
% Simulates response of GIF model to stimulus Stim
%
%  Inputs:  gifprs - param structure with fields: filt, vleak, 
%                      tau, sig, iht, ih 
%           Stim - stimulus matrix (time goes along vertical axis)
%           dtsim - bin size (i.e. fraction of a frame to use for sampling
%                             model voltage) 
%
%  Dynamics: Filters stimulus with 'filt' to obtain input current,
%            upsamples this to fine binning (size dtsim), adds Gaussian
%            noise (stdev 'sig') to obtain an OH process.  Spikes occur
%            whenever voltage exceeds 1 and resets to 0 in the next time
%            bin. Post-spike current 'ih' added to injected current
%            whenever a spike occurs.
%
%  Outputs:  tsp - vector of spike times (in units of frames).
%            Vmem - voltage trace
%            Ispk - total spiking-dependent current (ie. sum of ih's)
%            Istm - total stimulus-driven current
%               (total current = Istm + Ispk + Inoise)
%
%  Dependency: rundynam_gif_mex.c


% -- Process Inputs --------------------------------------
if 1./dtsim ~= round(1/dtsim)
    error('dtsim (binsize) must evenly divide 1');
end
nbns = 1/dtsim;  % number of output bins per input bin
[slen,swid] = size(Stim); % stimulus size
rlen = slen*nbns;  % number of output bins

tau = gifprs.tau;      % time const of exponential leak
vleak = gifprs.vleak;  % leak current reversal potential
sig = gifprs.sig;      % stdev of injected current noise
vreset = 0;  % reset potential after spike;  no loss of generality 
vthr = 1;    % taking vreset=0 and vthr=1; (but must have vthr>vreset)

% --  Compute filtered resp to signal ----------------------
filt = gifprs.filt;  % stim kernel
[nk1,nk2] = size(filt);
if (swid == nk2)
    Iinj = sameconv(Stim,filt);  % convolve filter with stim
else
    error('Mismatch between stimulus and filter size');
end

% -- Compute finely-sampled post-spike (h) current -----------------
ih = gifprs.ih;  % post-spike current
iht = gifprs.iht; % time sampling of post-spike current
ihthi = [dtsim:dtsim:max(iht)]';  % time points for sampling
ihhi = interp1(iht, ih, ihthi, 'linear', 0)*dtsim;

% ----- Run dynamics  ---------
decay1 = exp(-dtsim/tau);   % decay of V during 1 time bin
decay2 = (tau/dtsim)*(1-exp(-dtsim/tau)); % decay of Iinj contribution 
                                          % to V during 1 time bin
rndseed = round(rand*2.1e9);  % seed for random num generator
[Sp, Vmem, Ispk] = rundynam_gif_mex(Iinj, ihhi, vleak, vthr, vreset, ...
   sig, decay1, decay2, nbns, rndseed);  % mex code for running dynamics
tsp = find(Sp)*dtsim;
Ispk = Ispk/dtsim;

if nargout > 3 % Compute stim-dep current, if desired
    Istm = fastinterp(Iinj,100);  % Stim-based current only
end

%===============================================================
% Matlab equivalent, for error checking 
if 0
    [Sp2, Vmem2, Ispk2,Itot] = rundynam_gif(Iinj, ihhi, vleak, ...
        vthr, vreset, sig, decay1, decay2, nbns, rndseed);
    % note: Itot = Istm + Ispk + Inse
    Ispk2 = Ispk2/dtsim;
    Istm = fastinterp(Iinj,100);  % Stim-based current only
    plot([Itot [Istm+Ispk2]]);  % Should match if sig = 0;
end
% ==============================================================
