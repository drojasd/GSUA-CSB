function y = gsua_ua(M,T,varargin)
% Function for uncertainty analysis
%
% Y=gsua_ua(M,T)
% Parameters:
% M <-- design matrix from gsua_dmatrix function
% T <-- summary table from gsua_dataprep function
% Output:
% Y <-- result of Monte-Carlo simulation
% Additional features:
% The function automatically apply a Monte-Carlo filtering over results and
% present it as figures (gsua_MCF function). It is possible to perform UA
% over multiple model outputs. To obtain extrapolated results, use the
% paired feature 'xdata'. To compare results with a specific output use
% the paired feature 'ynom' (ynom length must coindice with xdata length).
% To avoid parallel computing (no speed up), use the paired feature 
% 'parallel',false.
% Y=gsua_ua(M,T,'xdata',xdata,'ynom',ynom,'parallel',false)
p=inputParser;
defNum=[];

addRequired(p,'M',@isnumeric);
addRequired(p,'T',@istable);
addParameter(p,'xdata',defNum,@isnumeric);
addParameter(p,'ynom',defNum,@isnumeric);
addParameter(p,'parallel',true,@islogical);

parse(p,M,T,varargin{:})
M=p.Results.M;
T=p.Results.T;
xdata=p.Results.xdata;
ynom=p.Results.ynom;
parallel=p.Results.parallel;

try
    TP=T.Properties.CustomProperties;
catch 
    TP=load('ATable.mat');
    TP=TP.Table2;
end
kind=TP.Kind;
out=TP.output;
nout=size(out,2);
vars=TP.Vars;


if strcmp('mat',kind)
    domain=TP.Domain;
    if isempty(xdata)
        xdata=linspace(domain(1),domain(2),1+domain(2)-domain(1));
    end
    if isempty(ynom)
        ynom=gsua_deval(T.Nominal',T,xdata);
    end
else
    if isempty(xdata)||isempty(ynom)
        disp('Generating xdata and ynom values')
        ynom=gsua_deval(T.Nominal',T,[]);
        xdata=0:size(ynom,2)-1;
    end
end

y=gsua_pardeval(M,T,xdata,parallel,[0,size(M,1)]);

figure(1)
clf
gsua_plot('UncertaintyAnalysis',T,y,xdata,ynom)
for i=1:nout
    figure('name',vars{out(i)});
    clf
    [~]=gsua_MCF(T,M,y(:,:,i),ynom(i,:))
end
end