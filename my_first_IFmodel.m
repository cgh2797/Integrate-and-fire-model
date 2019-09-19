%IF model
r=100;%Mohm
c= 200;%pF
e=-70;%mV
v=-60;%mV
vreset=e;
dt=0.01;%msec
v=e;%v0=e
i=2000;%?? ÀÌ°Å¹¹Áö
v_trace=[];
dur=dt*1000;
vthresh=-55;
niter=1000;
vm = zeros(1,niter);
vm(1)=v;
for idx= 2: niter
dvdt=(e-v+r*i)/(r*c) ;%difference eq;
v=v+dt.*dvdt;
     if v > vthresh
         v=vreset;
     end
     vm(idx)=v;
end
plot(vm)
title('my first IF-model')