function [Y,t] = sens_montecarlo(Par,M,xdata,parallel,init)
% [Y,t] = sens_montecarlo(model,M,init)
% Uncertainty analysis of a dynamical system (Simulink model) by generation of N sets of
% parameters and N simulations using the Monte Carlo method.
%
% model     Simulink model name. Example: 'Pendulum'. The Simulink model has to be set correctly:
%           (1) the name of parameters are p(1), p(2),...; (2) connect an "Out block" to the output;
%           (3) in "Configuration parameters | Data Import/Export" check these options: Time (tout),
%           Output (yout), Save simulation output as single object; (4) use fixed-step
%           solver; (5) if there are problems with rate transition, then fix the problem in
%           Configuration | Diagnostics | Sample time.
% M         Sample matrix (NxNp) of parameters (one parameter by column and one sample by row)
% Y         Matrix (NxNd) with time responses in rows to every of N set of parameters
% t         Time row-vector: t(i),i=1,2,...,Nd
% init      (Optional) Information for compute of percent progress and remaining time: [i0, Nsim], where i0 is the 
%           initial index of loop and Nsim is the total number of Montecarlo simulations
%
% (c) Carlos Mario Vélez S. 2015
% Universidad EAFIT, http://www.eafit.edu.co
% Medellín, Antioquia, Colombia
% E-mail: cmvelez@eafit.edu.co
% https://plus.google.com/+CarlosMVelezS/about
sim_progress='off';
try
    AT=Par.Properties.CustomProperties;
catch
    AT=load('ATable.mat');
    AT=AT.Table2;
end
model=AT.Kind;
output=AT.output;
nout=size(output,2);
Np=size(M,2);
load_system(model)

global t0
if nargin<5
    rem_compute=0;
else
    rem_compute=1;
    if init(1)==0
        t01=clock; t0=t01(4)*3600+t01(5)*60+t01(6);
    end
end


initT= str2double(get_param(model, 'StartTime'));
fixedStep = str2double(get_param(model, 'FixedStep'));
stopTime = str2double(get_param(model, 'StopTime'));
t=initT:fixedStep:stopTime;
if isempty(xdata)
    xdata=initT:stopTime;
end
Nd = length(xdata);
N = size(M,1); % Sample size
% Simulation with every set of parameters (changing all parameters) and storage in matrix Y (a column for every parameter set)
Y = zeros(N,Nd,nout);
out(N) = Simulink.SimulationOutput;
in(N) = Simulink.SimulationInput;
for i=1:N
    in(i)=in(i).setModelName(model);
    for w=1:Np
        in(i)=in(i).setVariable(Par.Properties.RowNames{w},M(i,w),'Workspace',model);
    end

end
if parallel
try
    
    out = parsim(in, 'ShowSimulationManager', sim_progress,'ShowProgress','off','UseFastRestart','off');
catch
    switch sim_progress
        case 'on'
        disp('SimulationManager has not been detected');
    end
    out = parsim(in,'ShowProgress',sim_progress,'UseFastRestart','off');
end
else
    for i=1:N
        out(i) = sim(in(i));
    end
end
failed=0;
% if N>1
    for i = 1:N
        try
            ytemp = interp1(out(i).get('tout'),out(i).get('yout'),xdata')'; % Y is a matrix of NxNd
            for j=output
                Y(i,:,j)=ytemp(j,:);
            end
        catch
            failed=failed+1;
        end
    end
% else
%     try
%         Y = interp1(out(1).get('tout'),out(1).get('yout'),xdata)'; % Y is a matrix of NxNd
%     catch
%         failed=failed+1;
%     end
% end
if failed>1
    disp({failed,'simulations have failed'})
end

if rem_compute==1 && N>1
        rem_percent=max(0,(1-(N+init(1))/init(2))*100);
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

