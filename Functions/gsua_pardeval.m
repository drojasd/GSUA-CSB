function y = gsua_pardeval(M,Table2,xdata,parallel,init)
%GSUA_PARDEVAL Evaluate the model behind a GSUA table for a batch of parameter sets
%
%   y = gsua_pardeval(M,T,xdata,parallel)
%   y = gsua_pardeval(M,T,xdata,parallel,init)
%
%   Batched counterpart of GSUA_DEVAL: evaluates the model behind T once
%   per row of the design/sample matrix M (optionally with parfor), and
%   assembles the results into an output array. Used internally by
%   GSUA_SA, GSUA_UA and GSUA_CSB to run the many simulations required
%   for sensitivity/uncertainty analysis.
%
%   Inputs:
%     M        <-- N x Np design matrix, one parameter set per row
%                  (from gsua_dmatrix)
%     T        <-- summary table from gsua_dataprep (or its
%                  CustomProperties struct directly)
%     xdata    <-- domain at which to evaluate/interpolate the model
%                  output
%     parallel <-- (logical) evaluate rows of M with parfor
%     init     <-- (optional) [i0,Nsim] progress-tracking info passed to
%                  gsua_timer-style reporting: i0 is the starting index
%                  of this call within a larger run, Nsim is the total
%                  number of simulations across that run
%
%   Outputs:
%     y <-- N x numel(xdata) (or N x numel(xdata) x nout for
%           multi-output models) array of model outputs
%
%   See also GSUA_DEVAL, GSUA_SA, GSUA_UA, GSUA_CSB, GSUA_DMATRIX.
try
    Table=Table2.Properties.CustomProperties;
catch
    Table=load('ATable.mat');
    Table=Table.Table2;
end
kind=Table.Kind;
if kind==1
    [y,~]=sens_montecarlo(Table2,M,xdata,parallel,init);
    return
end

global t0
if nargin<5
    rem_compute=0;
else
    rem_compute=1;
    if init(1)==0
        t01=clock; t0=t01(4)*3600+t01(5)*60+t01(6);
    end
end

ni=Table.NumVars;
fun=Table.Solver; 
out=Table.output;
N = size(M,1);
nout=size(out,2);
if parallel
    parforArg=Inf;
else
    parforArg=0;
end

switch kind
    case 3
        if nout==1 
            y=zeros(size(M,1),size(xdata,2));

            parfor (i=1:size(M,1),parforArg)
                try
                    y(i,:)=deval(fun(M(i,1:ni),M(i,ni+1:end)),xdata,out);
                catch
                end
            end
        else
            y=zeros(size(M,1),size(xdata,2),nout);
            parfor (i=1:size(M,1),parforArg)
                try
                    y_temp2=deval(fun(M(i,1:ni),M(i,ni+1:end)),xdata,out);
                    for j=1:nout
                        y(i,:,j)=y_temp2(j,:);
                    end
                catch
                end
            end
        end
    case 2
        y=fun(M);       
    case 4
        y=zeros(size(M,1),size(xdata,2));
        if ~parallel
            parfor i=1:size(M,1)
                y(i,:)=gsua_intrp(fun(M(i,1:ni),M(i,ni+1:end)),xdata,out);
            end
        else
            siz=size(M,1);
            div=divisors(siz);
            mi=min(div(div>=20));
            step=siz/mi;
            for i=1:mi
                sec=(i-1)*step+1:i*step;
                resp=fun(M(sec,1:ni),M(sec,ni+1:end));
                y(sec,:)=gsua_intrp(resp,xdata,out);
            end
        end
    case 5
        domain=Table.Domain;
        opt=Table.copt;
        if nout==1 
            y=zeros(size(M,1),size(xdata,2));

            parfor (i=1:size(M,1),parforArg)
                try
                    y(i,:)=gsua_intrp(fun(M(i,:),domain,opt),xdata,out);
                catch
                end
            end
        else
            y=zeros(size(M,1),size(xdata,2),nout);
            parfor (i=1:size(M,1),parforArg)
                try
                    y(i,:,:)=gsua_intrp(fun(M(i,:),domain,opt),xdata,out);
                catch
                end
            end
        end
    case 6
        y=zeros(size(M,1),1);
        parfor (i=1:size(M,1),parforArg)
            y(i)=fun(M(i,:));
        end
            
end   
if rem_compute==1
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
