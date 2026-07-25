function [cost_data,cost_band,P5,P50,P95] = gsua_covmetric(Y,ydata,varargin)
%GSUA_COVMETRIC Percentile-band coverage/tightness metrics for uncertainty-analysis results
%
%   [cost_data,cost_band,P5,P50,P95] = gsua_covmetric(Y,ydata)
%   [cost_data,cost_band,P5,P50,P95] = gsua_covmetric(Y,ydata,'margin',m,'alpha',a,'plevels',[lo hi])
%
%   Reduces a Monte-Carlo uncertainty-analysis output Y (from GSUA_UA / GSUA_PARDEVAL) to a
%   low/median/high percentile band per output signal and time step, then scores that band with
%   two costs built on the same margin-normalized, correlation-penalized MSE used throughout the
%   toolbox (see GSUA_COSTF): a value below 1 means "within the acceptable margin", a value at or
%   above 1 means "not yet within tolerance". This gives the semi-automation identification
%   routine a single, directly interpretable stopping signal instead of ad hoc envelope checks.
%
%   cost_data (accuracy)  <-- distance from the band's median curve to ydata. High means the
%       identified region's central tendency does not actually track the observed data.
%   cost_band (precision) <-- distance between the low and high percentile curves themselves,
%       independent of ydata. Low means repeated estimations/samples landed in a tight region,
%       i.e. the identification has converged, REGARDLESS of how well that region fits the data.
%       This distinction matters: a model can have a genuinely thin, converged identifiable
%       region that still does not perfectly track noisy or structurally-imperfect data. That is
%       a model-fit limitation, not an identifiability problem -- a semi-automation routine
%       should recognize low cost_band as "stop refining, this has converged" even when
%       cost_data never reaches 0.
%
%   Inputs:
%     Y       <-- reps x len (x nout) array of simulated trajectories, one row per Monte-Carlo
%                 run/estimation (output of gsua_ua/gsua_pardeval). Runs containing any
%                 non-finite value are excluded from the percentile band, with a warning.
%     ydata   <-- nout x len (or 1 x len for a single output) reference/experimental data
%     margin  <-- (optional) relative tolerance used to build the shared normalization
%                 regulator, same dual-mode convention as gsua_rcostf/gsua_costfMulti (values
%                 <1 are treated as a fraction and converted to 1+margin). Default: 0.1
%     alpha   <-- (optional) exponent applied to each per-signal cost term, same convention as
%                 gsua_costf. Default: 2
%     plevels <-- (optional) [low,high] percentile levels defining the band. Default: [5 95]
%
%   Outputs:
%     cost_data <-- scalar normalized cost of the band median vs ydata
%     cost_band <-- scalar normalized cost of the high percentile curve vs the low percentile
%                   curve (band tightness)
%     P5,P50,P95 <-- nout x len low/median/high percentile curves actually used (named for the
%                    default plevels but reflecting whatever plevels was passed; returned for
%                    plotting/reuse, e.g. with gsua_plot)
%
%   Example (semi-automation routine, Phase 1 reachability / Phase 4 stopping check):
%     Y = gsua_ua(M,T,'xdata',xdata,'ynom',ydata,'parallel',true);
%     [cost_data,cost_band] = gsua_covmetric(Y,ydata,'margin',0.1);
%     reachable = cost_data < 1;               % Phase 1
%     converged = cost_data < 1 && cost_band < 1;  % Phase 4
%
%   See also GSUA_COSTF, GSUA_UA, GSUA_IA, GSUA_MEDIANCI, GSUA_PE.

p=inputParser;
addRequired(p,'Y',@isnumeric);
addRequired(p,'ydata',@isnumeric);
addParameter(p,'margin',0.1,@isnumeric);
addParameter(p,'alpha',2,@isnumeric);
addParameter(p,'plevels',[5 95],@(x) isnumeric(x) && numel(x)==2 && x(1)<x(2));

parse(p,Y,ydata,varargin{:})
Y=p.Results.Y;
ydata=p.Results.ydata;
margin=p.Results.margin;
alpha=p.Results.alpha;
plevels=p.Results.plevels;

if ismatrix(Y)
    Y=reshape(Y,size(Y,1),size(Y,2),1);
end
if isvector(ydata)
    ydata=reshape(ydata,1,[]);
end

nout=size(Y,3);
len=size(Y,2);

valid=all(all(isfinite(Y),2),3);
nBad=sum(~valid);
if nBad>0
    warning('gsua_covmetric:invalidRuns',...
        '%u of %u runs contain non-finite values and were excluded from the percentile band',...
        nBad,size(Y,1))
end
Y=Y(valid,:,:);
if size(Y,1)<2
    error('gsua_covmetric:tooFewRuns','Not enough finite runs left to compute a percentile band')
end

P5=zeros(nout,len); P50=zeros(nout,len); P95=zeros(nout,len);
for j=1:nout
    P5(j,:)=prctile(Y(:,:,j),plevels(1),1);
    P50(j,:)=prctile(Y(:,:,j),50,1);
    P95(j,:)=prctile(Y(:,:,j),plevels(2),1);
end

cost_data=local_bandcost(ydata,P50,len,margin,alpha);
cost_band=local_bandcost(P5,P95,len,margin,alpha);

end

function cost=local_bandcost(a,b,len,margin,alpha)
%Correlation-penalized, margin-normalized MSE cost, same shape as GSUA_COSTF, with the
%isnan(corr2(...)) guard already used in GSUA_COSTFMULTI so a perfectly flat reference/band
%(the converged, thin-band case this function exists to detect) does not silently propagate
%NaN through the metric.
margin=abs(margin);
if margin<1
    margin=margin+1;
end
inputs=size(a,1);
regulator=sum((a-a*margin).^2,2)/len;
regulator(regulator==0)=eps;
cost=sum((a-b).^2,2)/len./regulator;
for i=1:inputs
    r=corr2(a(i,:),b(i,:));
    if isnan(r)
        r=1;
    end
    cost(i)=((2-r)*cost(i))^alpha;
end
cost=sum(cost)/inputs;
end
