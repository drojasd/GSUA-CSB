function t1 = gsua_timer(rem_compute,init)
%GSUA_TIMER Print elapsed/remaining-time progress for long simulation loops
%
%   gsua_timer(rem_compute,init)
%   t1 = gsua_timer(rem_compute,init)
%
%   Internal progress-reporting helper used by GSUA_PE and
%   GSUA_LIKELIHOOD. Call once with rem_compute=0 to start the shared
%   timer (stored in the global t0), then repeatedly with rem_compute=1
%   to print progress percentage, elapsed time, and estimated remaining
%   time to the command window.
%
%   Inputs:
%     rem_compute <-- 0 to (re)start the timer, 1 to print a progress
%                     update
%     init        <-- [i,N] current iteration index and total number of
%                     iterations
%
%   Outputs:
%     t1 <-- current clock reading in seconds-of-day (only set when
%            rem_compute is 1)
%
%   See also GSUA_PE, GSUA_LIKELIHOOD.

global t0
if init(1)==0
    t01=clock; t0=t01(4)*3600+t01(5)*60+t01(6);
end
if rem_compute==1
    rem_percent=max(0,(1-init(1)/init(2))*100);
    t11=clock; t1=t11(4)*3600+t11(5)*60+t11(6);
    dt=t1-t0;
    dt_h=floor(dt/3600); dt_m=floor((dt-dt_h*3600)/60);dt_s=dt-dt_h*3600-dt_m*60;
    dte=dt/(1-rem_percent/100);
    dte_h=floor(dte/3600); dte_m=floor((dte-dte_h*3600)/60);dte_s=dte-dte_h*3600-dte_m*60;
    te=t0+dte;
    te_h=floor(te/3600); te_m=floor((te-te_h*3600)/60); te_s=te-te_h*3600-te_m*60;
    tr=max(0,te-t1);
    tr_h=floor(tr/3600); tr_m=floor((tr-tr_h*3600)/60); tr_s=tr-tr_h*3600-tr_m*60;
    clc
    disp(['Progress: ' num2str(100-floor(rem_percent)) '%'])
    disp(['Estimated processing time (h:m:s): ' num2str(dte_h) ':' num2str(dte_m) ':' num2str(floor(dte_s))])
    disp(['Remaining time (h:m:s): ' num2str(tr_h) ':' num2str(tr_m) ':' num2str(floor(tr_s))])
    disp(['Elapsed time (h:m:s): ' num2str(dt_h) ':' num2str(dt_m) ':' num2str(floor(dt_s))])      
    disp(['Estimated stop time (h:m:s): ' num2str(te_h) ':' num2str(te_m) ':' num2str(floor(te_s))])
    disp(['Number of simulations: ' num2str(init(2))])
end

end