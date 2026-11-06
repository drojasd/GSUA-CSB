function [Table,res,IP] = gsua_pe(Table,xdata,ydata,varargin)
% Parameter estimation function
% [T,res] = gsua_pe(T,xdata,ydata,N)
% Parameters:
% T     <-- summary table from gsua_dataprep function
% xdata <-- array of points where the model will be evaluated
% ydata <-- array with expected model output
% Outputs:
% T   <-- summary table with parameter estimation results. Besides T.Est<solver> (Np x N, columns
%         sorted best-first), it now carries T.IP<solver> (Np x N, the initial point actually used
%         for each run, same column order) and the CustomProperty Retries<solver> (1 x N, how many
%         silent restarts from a random start point each estimate cost), so the retry behaviour is
%         traceable instead of silent.
% res <-- cost functions for each estimation
% Also, for the gradient solvers ('fmincon','lsqc','lsqn'), if you do not pass an 'opt' with a
% TypicalX set, gsua_pe fills TypicalX with the nominal factor magnitudes (max(abs(T.Nominal),eps)),
% which greatly improves finite-difference scaling when factors span many orders of magnitude.
% Additional paired features:
% 'N',N     <-- number of parameter estimations
% 'Multistart',k <-- activate multistart feature for lsqcurvefit optimizer,
% perform k parameter estimations in each cycle
% 'solver',{'lsqc','lsqn','ga','particle','psearch','surrogate','annealing'}
% <-- allows to choose among several matlab optimizers. Default: lsqc
% 'opt',optimoptions(optimizer,...) <-- allows to configure the respective
% matlab optimizer, optimizer must match with the real matlab optimizer
% name.
% 'ipoint',point <-- allows to give to the optimizer the initial points for
% estimations. point must have the same number of columns as factors to
% estimate and the same number of rows as N.
% 'Show',{'off','on'}. If activated, the function plot optimizer results
% Default: off.
% 'save',{false,true} <-- if true, writes Estimations.mat to the CURRENT folder after
% each estimation (a crash-resume checkpoint), OVERWRITING any existing
% Estimations.mat there. Default: false (it was true before v1.15).
% 'margin',m <-- relative tolerance band for the correlation-penalized cost (GSUA_COSTF), used by
% the 'ga','particle','psearch','surrogate','annealing' and 'fmincon' solvers. The value you pass
% is a RELATIVE band on the residual: with a good-shape fit (correlation ~ 1) the cost drops below
% 1 exactly when the residual RMS is within a fraction m of the data's own RMS -- so 'margin',0.3
% means "accept, as cost < 1, fits whose residual is within ~30% of the data scale". Internally the
% regulator is built with m+1 (so 0.3 is used as 1.3); you ALWAYS pass the raw band m, and
% GSUA_NOISEFLOOR/GSUA_BOOTSTRAP take the same raw m. The regulator does not depend on the
% parameters, so margin only RESCALES the cost: it does NOT move the optimum, only the
% interpretation of cost thresholds (cost < 1, the noise floor, the fit-quality cutoff). m = 0
% selects plain MSE; m < 0 selects the Gaussian negative-log-likelihood path (GSUA_LIKELIHOOD).
% Default: 0.
% 'alpha',a <-- exponent applied to each per-output term in GSUA_COSTF. Default: 2.
% [T,res] = gsua_pe(T,xdata,ydata,'Show','on','solver','particle','N',3)
% [T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.3)  % cost<1 ~ residual within 30%
%
% See also GSUA_DMATRIX, GSUA_DEVAL, GSUA_COSTF, GSUA_IA, GSUA_DIA.

p=inputParser;

dfN=1;
dfMS=1;
dfIP=[];
dfOS='lsqc';
dfOpt=[];

validOS={'lsqc','lsqn','ga','particle','psearch','surrogate','annealing','fmincon'};
checkOS = @(x) any(validatestring(x,validOS));
checkIP = @(x) (isnumeric(x) || isempty(x));
defaultShow='off';
validShow={'off' 'on' 'on2'};
checkShow = @(x) any(validatestring(x,validShow));
dfmarg=0;
dfalpha=2;
dfSave=false;   % opt-in since v1.15: writing Estimations.mat silently overwrote user files
dfTimer=true;

addRequired(p,'Table',@istable);
addRequired(p,'xdata',@isnumeric);
addRequired(p,'ydata',@isnumeric);
addParameter(p,'N',dfN,@isnumeric);
addParameter(p,'Multistart',dfMS,@isnumeric);
addParameter(p,'ipoint',dfIP,checkIP);
addParameter(p,'solver',dfOS,checkOS);
addParameter(p,'opt',dfOpt);
addParameter(p,'Show',defaultShow,checkShow);
addParameter(p,'A',dfOpt);
addParameter(p,'B',dfOpt);
addParameter(p,'Aeq',dfOpt);
addParameter(p,'Beq',dfOpt);
addParameter(p,'nonlcon',dfOpt);
addParameter(p,'margin',dfmarg,@isnumeric);
addParameter(p,'alpha',dfalpha,@isnumeric);
addParameter(p,'save',dfSave,@islogical);
addParameter(p,'timer',dfTimer,@islogical);

parse(p,Table,xdata,ydata,varargin{:})
T=p.Results.Table;
xdata=p.Results.xdata;
ydata=p.Results.ydata;
N=p.Results.N;
MS=p.Results.Multistart;
IP=p.Results.ipoint;
OS=p.Results.solver;
Opt=p.Results.opt;
show=p.Results.Show;
A=p.Results.A;
B=p.Results.B;
Aeq=p.Results.Aeq;
Beq=p.Results.Beq;
nonlcon=p.Results.nonlcon;
margin=p.Results.margin+1;
saver=p.Results.save;
costAlpha = p.Results.alpha;
timer=p.Results.timer;

% Record the exact margin/alpha this call will score with, so a downstream noise-floor
% calibration (GSUA_NOISEFLOOR) can recompute the identical regulator/cost formula without the
% caller having to re-specify (and risk mismatching) them. Margin is stored already offset by +1,
% matching the variable actually used in the regulator formula below -- only meaningful, and only
% set, for the branch that actually scores with GSUA_COSTF.
if any(strcmp(OS,{'ga','particle','psearch','surrogate','annealing','fmincon'})) && margin ~= 1 && margin > 0
    try
        Table = addprop(Table,{'Margin','Alpha'},{'table','table'});
    catch
    end
    try
        Table.Properties.CustomProperties.Margin = margin;
        Table.Properties.CustomProperties.Alpha = costAlpha;
    catch
    end
end

lb=T.Range(:,1);
ub=T.Range(:,2);
Np=size(T,1);

err_limit=100;

if MS>1
    OS='Multistart';
    if timer, disp('only lsqc supports Multistart ... Running lsq optimization'), end
end
if (isempty(IP) || size(IP,1)~=N) && any(strcmp(OS,{'lsqc','Multistart','lsqn','fmincon','psearch','annealing'}))
    if timer, disp('Generating a valid matrix for estimations'), end
    [IP,~]=gsua_dmatrix(T,N);
end

% Default TypicalX for the gradient solvers when the user did not set one. With factors spanning
% many orders of magnitude (e.g. 1e-4 to 1e4) a single global finite-difference step behaves poorly;
% TypicalX = the nominal magnitudes rescales the differences. Only filled in when absent, so an
% explicit user optimoptions is never overridden.
gradSolver = '';
switch OS
    case 'fmincon', gradSolver = 'fmincon';
    case 'lsqc',    gradSolver = 'lsqcurvefit';
    case 'lsqn',    gradSolver = 'lsqnonlin';
end
if ~isempty(gradSolver)
    if ismember('Nominal', T.Properties.VariableNames)
        typX = max(abs(T.Nominal(:).'),eps);
    else
        typX = max(abs(mean(T.Range,2).'),eps);
    end
    if isempty(Opt)
        Opt = optimoptions(gradSolver);
    end
    try
        if isempty(Opt.TypicalX)
            Opt = optimoptions(Opt,'TypicalX',typX);
        end
    catch
        % options object without a settable TypicalX: leave Opt untouched
    end
end

x=zeros(N,Np);
res=inf(N,1);
retries=zeros(N,1);   % number of silent restarts (random IP) spent on each estimate
warnedZeroReg=false;  % so the zero-cost-scale warning (below) fires at most once
lastME=[];   % remember the last genuine estimation error, for the diagnostics below

% Pre-flight: the most common reason gsua_pe "returns only Inf" is a model that does
% not evaluate to finite output (bad model, wrong sol orientation, xdata outside the
% integration domain, ...). gsua_deval swallows the model error and returns Inf, so the
% optimizer "succeeds" with an Inf cost and the user is left with an all-Inf res and no
% cause. Catch it here, before spending the optimizer budget, with a message that says
% what to inspect.
if ismember('Nominal', T.Properties.VariableNames)
    lastwarn('');
    y_nom_check = gsua_deval(T.Nominal', T, xdata);
    [wmsg, ~] = lastwarn;
    if ~all(isfinite(y_nom_check(:)))
        hint = 'Check that your model returns finite output over xdata; inspect it with  gsua_deval(T.Nominal'',T,xdata).';
        if isempty(wmsg)
            error('GSUA:gsua_pe:ModelNotFinite', ...
                'The model returns non-finite output at the nominal parameters, so no estimation can succeed.\n%s', hint);
        else
            error('GSUA:gsua_pe:ModelNotFinite', ...
                'The model returns non-finite output at the nominal parameters, so no estimation can succeed.\n%s\nLast model warning: %s', hint, wmsg);
        end
    end
end
% try
%     TP=T.Properties.CustomProperties;   
% catch
%     TP=load('ATable.mat');
%     TP=TP.Table2;
% end
% if any(strcmp(OS,{'lsqc','lsqn'})) && N>1
%     
%     parfor i=1:N
%         try
%             switch OS
%                 case 'lsqc'
%                      [x(i,:),res(i)] = lsqcurvefit(@(pars,xaxis) gsua_deval(pars,TP,xaxis),IP(i,:),xdata,ydata,lb,ub,Opt);
%                 case 'lsqn'
%                     obj=@(pars) gsua_deval(pars,TP,xdata)-ydata;  
%                     [x(i,:),res(i)] = lsqnonlin(obj,IP(i,:),lb,ub,Opt);
%             end
%         catch ME
%             warning(ME.message)  
%         
%         end
%     end    
% else
if timer
    gsua_timer(0,[0,N]);
end
    for i=1:N
        err_counter=0;
        if timer, disp(['Estimation ' num2str(i)]), end
        done=true;
        attempt=1;
    while done
    try
    switch OS
        case 'Multistart'
            problem = createOptimProblem(OS,'x0',IP(i,:),'objective',@(pars,xaxis) gsua_deval(pars,T,xaxis),...
            'lb',lb,'ub',ub,'xdata',xdata,'ydata',ydata);
            ms = MultiStart(UseParallel',1);
            [x(i,:),res(i)] = run(ms,problem,multi);           
        case 'lsqc'
            [x(i,:),res(i)] = lsqcurvefit(@(pars,xaxis) gsua_deval(pars,T,xaxis),IP(i,:),xdata,ydata,lb,ub,Opt);
        case 'lsqn' 
             obj=@(pars) (gsua_deval(pars,T,xdata)-ydata);   
             [x(i,:),res(i)] = lsqnonlin(obj,IP(i,:),lb,ub,Opt);
        case {'ga','particle','psearch','surrogate','annealing','fmincon'}
             if margin ~= 1 && margin>0
                 inputs=size(ydata,1);
                 len=length(xdata);
                 % omitnan + per-row non-NaN count (not the fixed len) so a row with
                 % gappy real-world data isn't silently deflated relative to a
                 % fully-populated row, and a NaN doesn't propagate into regulator
                 % (which would poison gsua_costf's cost for every row, not just the
                 % gappy one).
                 regulator=sum((ydata-ydata*margin).^2,2,'omitnan')./sum(~isnan(ydata),2);
                 % Guard a zero/undefined cost scale (an all-zero or constant data row -- e.g. a 0-1
                 % fraction row whose synthetic draw came out all zeros): gsua_costf would divide by
                 % it and hand fmincon an undefined objective, failing every start without saying
                 % which row. Set it to Inf so that row contributes 0, and name the row(s) once.
                 badReg = ~(regulator>0 & isfinite(regulator));
                 if any(badReg)
                     if ~warnedZeroReg
                         warning('GSUA:gsua_pe:ZeroRegulator', ...
                             ['Output row(s) %s have a zero/undefined cost scale (all-zero or ' ...
                              'constant data); excluding them from the cost. For a joint count+' ...
                              'fraction fit, set a per-row noise model in gsua_bootstrap.'], ...
                             mat2str(find(badReg(:).')));
                         warnedZeroReg = true;
                     end
                     regulator(badReg) = Inf;
                 end
                 obj=@(pars) gsua_costf(inputs,regulator,len,ydata,gsua_deval(pars,T,xdata),costAlpha);
             else
                 if margin <0
                     desv = abs((ydata*margin).^2);
                     obj =@(pars) sum(log(2*pi*desv) +(ydata-gsua_deval(pars,T,xdata)).^2./desv,'all','omitnan')/2;
                 else
                    obj=@(pars) gsua_mse(ydata,gsua_deval(pars,T,xdata));
                 end
             end
            switch OS
                case 'ga'
                    [x(i,:),res(i)] = ga(obj,Np,A,B,Aeq,Beq,lb,ub,nonlcon,Opt);
                case 'particle'
                    [x(i,:),res(i)] = particleswarm(obj,Np,lb,ub,Opt);    
                case 'psearch'
                    [x(i,:),res(i)] = patternsearch(obj,IP(i,:),A,B,Aeq,Beq,lb,ub,nonlcon,Opt);
                case 'surrogate'
                    [x(i,:),res(i)] = surrogateopt(obj,lb,ub,Opt);
                case 'annealing'
                    [x(i,:),res(i)] = simulannealbnd(obj,IP(i,:),lb,ub,Opt);
                case 'fmincon'
                    [x(i,:),res(i)] = fmincon(obj,IP(i,:),A,B,Aeq,Beq,lb,ub,nonlcon,Opt);
            end
    end
    done=false;
    if saver
        save('Estimations.mat','x','res');
    end
    catch ME
        err_counter=err_counter+1;
        if err_counter>err_limit
            warning('GSUA:gsua_pe:ErrorLimit', ...
                ['Reached the %d-error limit on estimation %d of %d; %d estimation(s) left unattempted.\n' ...
                 'Last error (%s): %s'], err_limit, i, N, N-i, ME.identifier, ME.message);
            retries(i)=err_counter;
            [res,index] = sort(res);
            %res=res(res ~= inf);
            res=res';
            x=x(index(1:size(res,2)),:);
            Table.(['Est' char(OS)]) = x';
            [Table,IP] = local_storeStarts(Table,OS,IP,retries,index,N);
            if all(~isfinite(res))
                error('GSUA:gsua_pe:AllEstimationsFailed', ...
                    ['Every estimation failed before any converged (error limit hit on #%d).\n' ...
                     'Check that your model returns finite output over xdata; inspect it with  gsua_deval(T.Nominal'',T,xdata).\n' ...
                     'Last error (%s): %s'], i, ME.identifier, ME.message);
            end
            return
        end
        switch ME.identifier
            case 'optimlib:lsqcurvefit:YdataSizeMismatchFunVal'
            if attempt==1
                ydata=ydata';
                attempt=2;
            else
                throw(ME)
            end
            case 'MATLAB:warning:invalidMessageType'
                throw(ME)
            otherwise
                lastME=ME;
                warning('GSUA:gsua_pe:RetryEstimate', ...
                    'Estimation %d failed (%s): %s -- retrying from a new random start point.', ...
                    i, ME.identifier, ME.message);
                IP(i,:)=gsua_dmatrix(T,1);
        end
    end
    end
    retries(i)=err_counter;
    if timer
        gsua_timer(1,[i,N]);
    end
    end
% end


[res,index] = sort(res);
%res=res(res ~= 0);
res=res';
x=x(index(1:size(res,2)),:);
Table.(['Est' char(OS)]) = x';
[Table,IP] = local_storeStarts(Table,OS,IP,retries,index,N);
if all(~isfinite(res))
    % Every estimation returned a non-finite cost. Fail loudly with the likely cause
    % instead of handing back an all-Inf res the caller has to decode.
    if isempty(lastME)
        error('GSUA:gsua_pe:AllEstimationsFailed', ...
            ['All %d estimations returned non-finite cost; none converged.\n' ...
             'Check that your model returns finite output over xdata; inspect it with  gsua_deval(T.Nominal'',T,xdata).'], N);
    else
        error('GSUA:gsua_pe:AllEstimationsFailed', ...
            ['All %d estimations returned non-finite cost; none converged.\n' ...
             'Check that your model returns finite output over xdata; inspect it with  gsua_deval(T.Nominal'',T,xdata).\n' ...
             'Last error (%s): %s'], N, lastME.identifier, lastME.message);
    end
end
if strcmp(show,'on')||strcmp(show,'on2')
    figure(1);
    clf
    [~]=gsua_eval(x',Table,xdata,ydata);
    if strcmp(show,'on2')
        gsua_ia(Table,x')
    end
end
end

function e = gsua_mse(A,B)
%GSUA_MSE Mean squared error between two equally sized arrays.
%   e = gsua_mse(A,B) returns mean((A(:)-B(:)).^2), i.e. sum of squared
%   differences divided by the number of elements. This reproduces the
%   Image Processing Toolbox function IMMSE for double-valued data, so the
%   fmincon/ga/particleswarm plain-MSE objective in gsua_pe no longer
%   depends on that toolbox. A and B must be the same size; as with immse,
%   any NaN propagates to the result.
if ~isequal(size(A),size(B))
    error('GSUA:gsua_mse:SizeMismatch', ...
        'A and B must be the same size (got %s and %s).', ...
        mat2str(size(A)), mat2str(size(B)));
end
d = double(A) - double(B);
e = sum(d(:).^2) / numel(d);
end

function [Table,IP] = local_storeStarts(Table,OS,IP,retries,index,N)
%LOCAL_STORESTARTS Record, aligned to the SORTED estimates, the initial point actually used for each
%   run (Table.IP<solver>, Np x N, columns best-first like Est<solver>) and the number of silent
%   restarts it took (CustomProperty Retries<solver>, 1 x N), and return IP reordered to match
%   res/Est. This makes the previously silent "retry from a random point" traceable: you can see
%   which start produced each fit and how many attempts it cost.
if ~isempty(IP) && size(IP,1)==N
    IP = IP(index,:);
    Table.(['IP' char(OS)]) = IP.';
end
rprop = ['Retries' char(OS)];
try
    Table = addprop(Table,{rprop},{'table'});
catch
end
try
    Table.Properties.CustomProperties.(rprop) = reshape(retries(index),1,[]);
catch
end
end