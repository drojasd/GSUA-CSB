function [ParT,Y,J] = gsua_sa(M,ParT,varargin)
% Function for sensitivity analysis
%
% T=gsua_sa(M,T)
% Parameters:
% M <-- design matrix from gsua_dmatrix function
% T <-- summary table from gsua_dataprep function
% Output:
% T <-- summary table with sensitivity indices.
% Additional features:
% It is possible to request 2 extra positional outputs, in this order:
% Y <-- model output matrix (same as gsua_ua output), 2nd output
% J <-- vector of model outputs in scalar representation (MSE), 3rd output
% It is possible to choose between 7 SA methods using the paired feature
% 'SensMethod',method. Where method is one of the following:
% 'Saltelli'    <-- requires N*(Np/2+1) simulations
% 'Jansen'      <-- requires N*(Np/2+1) simulations
% 'Xiao'        <-- newest method (Default), requires N*(Np/2+1)
%                   simulations
% 'Sobol'       <-- classic method, requires N*(Np+1) simulations
% 'Shapley'     <-- Shapley effects (see SHAPLEY EFFECTS below); writes
%                   T.Sh / T.Sh_vec instead of T.Si / T.STi
% 'brute-force' <-- explores all possible combinations, N + 2Np*N^2
%                   simulations
% 'OAT'         <-- non-global method, requires N*(Np+1) simulations.
%                   Varies one factor at a time over its N sampled values
%                   with the others at T.Nominal; Si is each factor's share
%                   of the summed OAT variance of the cost J and STi that
%                   variance itself. No interactions, no Si_vec/STi_vec.
%
% WHAT EACH METHOD ESTIMATES
% 'Sobol' and 'Jansen' estimate the variance-based (Sobol') first-order
% and total indices of the output Y at each time (Si_vec, STi_vec).
% 'Xiao' (Xiao, Lu & Wang 2018, distance components / energy distance)
% estimates STi = E|YA-YABi|^a / E|YA-YB|^a with a = 'xAlpha' (default 1).
% With 'xAlpha',2 it coincides with the Jansen (variance-based) index; with
% the default 1 it is a moment-independent, distribution-based index, which
% is more robust but can rank the factors differently from Sobol'.
% 'Saltelli' scores the squared-deviation field (Y-ynom).^2 rather than Y
% itself, i.e. it measures which factors drive the squared departure from
% the nominal trajectory -- a different quantity from Sobol'/Jansen.
% For every method the scalar Si/STi are indices of the cost
% J = sum((Y-ynom).^2) (time and series summed).
%
% INDICES OUTSIDE [0,1]
% Sensitivity indices lie in [0,1] by definition, but their Monte Carlo
% ESTIMATES do not: with a small N an estimated STi can exceed 1 and an Si
% can be negative. That is sampling noise of the estimator, not an error,
% and it is deliberately not clipped (clipping would hide that N is too
% small). It is most frequent for 'Saltelli', occasional for 'Sobol' and
% 'Jansen', and rare for 'Xiao' with xAlpha=1. Increase N until
% sum(Si)/sum(abs(Si)) is stable near 1 (a ratio well below 1 means
% spurious negative Si, i.e. N is too small).
%
% SHAPLEY EFFECTS ('SensMethod','Shapley')
% Sh_i shares the output variance among the factors, splitting every
% interaction fairly between the factors involved (Owen 2014; Song, Nelson
% & Staum 2016). They always sum to 1 and, for independent factors,
% Si <= Sh_i <= STi. Unlike Sobol' indices they remain meaningful when the
% factors are DEPENDENT (correlated or structurally coupled, e.g.
% R0 = beta*S0/gamma). Choose the estimator with 'ShapleyType':
%   'Independent' (default) <-- factors sampled independently from T.Range
%       (as gsua_dmatrix does). Nested Monte Carlo over random permutations
%       (Song, Nelson & Staum 2016). Cost: N + m*(Np-1)*No*Ni simulations,
%       with N = size(M,1) used for Var(Y). Options: 'ShapleyPerm' m
%       (default N; when Np! <= m every permutation is used equally often),
%       'ShapleyNo' outer draws (default 1), 'ShapleyNi' inner draws
%       (default 3).
%   'Given' <-- a JOINT sample whose dependence is kept: the rows of M
%       (e.g. from gsua_dmatrix(Tia,N,'Method','Joint','JointType',
%       'SmoothBootstrap'); the default 'Bootstrap' only repeats pool
%       vectors) or 'Pool', an Np x nPool ensemble such as Tia.Est.
%       Nearest-neighbour estimator
%       of Broto, Bachoc & Depecker (2020); costs one simulation per
%       distinct sample vector (duplicates are dropped). Option
%       'ShapleyK' neighbours (default 3). It answers "how is the output
%       variance shared WITHIN this joint distribution": with an accepted-
%       estimate ensemble (Tia.Est) that is sensitivity inside the
%       identified region, not over the prior ranges. Use a few hundred
%       distinct vectors or more; small samples give noisy estimates.
% Like any Monte Carlo index, an estimated Sh_i can come out slightly
% negative at small N; the scalar T.Sh, built on the heavy-tailed cost J,
% needs a larger N than T.Sh_vec to settle.
% Results: T.Sh (Np x 1, indices of the cost J) and T.Sh_vec
% (Np x Nd x nSeries, indices of Y at each time). Plot them like the other
% indices, e.g. gsua_plot('Pie',T,T.Sh) or
% gsua_plot('FractionalSensitivityArea',T,T.Sh_vec,xdata).
%
% Also, it is possible to specify an output for indices estimation using
% the paired feature 'ynom'. If you want to run the analysis with a
% different resolution for the time-dependent sensitivity indices, you can
% provide the paired feature 'xdata', [init:step:final].
% Finally, you can avoid parallel speed up 
% for simulation process using the paired feature 'parallel', false.
% [ParT,Y,J] = gsua_sa(M,T,'SensMethod','Xiao','ynom',ynom)
% Note: ynom size must coincide with xdata size, where xdata is the next
% xdata=linspace(Domain(1),Domain(2),Domain(2)-Domain(1)+1) and Domain is
% the array of the model domain. If you are using the paired feature
% 'xdata', ignore the last advice.
%
% See also GSUA_UA, GSUA_DMATRIX, GSUA_PARDEVAL, GSUA_COSTFMULTI, GSUA_DEPTH.

p=inputParser;

defaultSensMethod='Xiao';
validSensMethod={'brute-force','Jansen','Sobol','OAT','Saltelli','Xiao','Shapley'};
checkSensMethod = @(x) any(validatestring(x,validSensMethod));
validShapleyType={'Independent','Given'};
checkPosInt = @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x);

defaultParallel=true;

defaultXdata=[];

defaultyexp=[];

checkxAlpha = @(x) (x>0 && x<=2);
dxAlpha=1;


addRequired(p,'M',@isnumeric);
addRequired(p,'ParT',@istable);
addParameter(p,'ynom',defaultyexp,@isnumeric);
addParameter(p,'SensMethod',defaultSensMethod,checkSensMethod);
addParameter(p,'parallel',defaultParallel,@islogical);
addParameter(p,'verbose',true,@islogical);
addParameter(p,'xdata',defaultXdata,@isnumeric);
addParameter(p,'xAlpha',dxAlpha,checkxAlpha);
addParameter(p,'bandwidth',dxAlpha,@isnumeric);
addParameter(p,'multiObjective',false,@islogical);
addParameter(p,'ShapleyType','Independent',@(x) any(validatestring(x,validShapleyType)));
addParameter(p,'ShapleyPerm',[],@(x) isempty(x) || checkPosInt(x));
addParameter(p,'ShapleyNo',1,checkPosInt);
addParameter(p,'ShapleyNi',3,@(x) checkPosInt(x) && x>=2);
addParameter(p,'ShapleyK',3,@(x) checkPosInt(x) && x>=2);
addParameter(p,'Pool',[],@isnumeric);

parse(p,M,ParT,varargin{:})

M=p.Results.M;
ParT=p.Results.ParT;
y_exp=p.Results.ynom;
% Canonical spelling, so 'xiao' or 'shapley' select the method instead of falling through.
SensMethod=validatestring(p.Results.SensMethod,validSensMethod);
shapType=validatestring(p.Results.ShapleyType,validShapleyType);
shapPerm=p.Results.ShapleyPerm;
shapNo=p.Results.ShapleyNo;
shapNi=p.Results.ShapleyNi;
shapK=p.Results.ShapleyK;
shapPool=p.Results.Pool;
Parallel=p.Results.parallel;
verbose=p.Results.verbose;
xdata=p.Results.xdata;
pod=p.Results.xAlpha;
margin=p.Results.bandwidth;


try
    PT=ParT.Properties.CustomProperties;
catch 
    PT=load('ATable.mat');
    PT=PT.Table2;
end

kind=PT.Kind;

% Multi-output models default to PER-SERIES time-dependent indices: the single-objective path below
% now carries a native series dimension, so Si_vec/STi_vec come back as (Np x Nd x nSeries). The
% former depth-based scalarization of several output series into ONE objective is opt-in via
% 'multiObjective',true (it returns scalar Si/STi only).
multiObjt = p.Results.multiObjective;

if isempty(xdata)
    if kind~=1
        domain=PT.Domain;
        if size(domain,2)==1
            xdata=domain;
            Parallel=false;
        else
            xdata=linspace(domain(1),domain(2),domain(2)-domain(1)+1);
        end
    else
        xdata=[];
    end
end

Par=ParT.Properties.RowNames;
fixed=PT.Fixed;



% if Parallel && exist('parpool','file')>0 && isempty(gcp('nocreate')) && (kind==1 || kind==3 || kind==5)
%     parpool local;
% elseif ~Parallel && exist('parpool','file')>0 && isempty(gcp('nocreate'))==0
%     delete(gcp('nocreate'))
% end

Np=size(Par,1);
Np2=size(Par,1)-sum(fixed);
N=size(M,1);
if mod(N,2)~=0
    N=N-1;
    warning('The number of combinations to evaluate must be an even number')
    disp("Rows in design matrix were reduced to: "+num2str(N))
    M=M(1:N,:);
end

if isempty(fixed)
    fixed=false(1,Np);
end

if ~multiObjt
    if isempty(y_exp)
        y_exp=gsua_deval(ParT.Nominal',ParT,xdata);
    end
    % Native series dimension: y_exp is (nSeries x Nd) and the batched Y is (N x Nd x nSeries) for a
    % multi-series (matrix) user model, or the usual 2-D (N x Nd) when nSeries==1. yexp3 is the
    % reference reshaped to (1 x Nd x nSeries) so it broadcasts against Y; summing a per-sample scalar
    % cost J over BOTH time and series keeps the scalar Si/STi defined while Si_vec/STi_vec stay
    % per-series (Np x Nd x nSeries). nSeries==1 reduces to the previous 2-D math exactly.
    nSeries = size(y_exp,1);
    Nd = size(y_exp,2);
    yexp3 = reshape(y_exp.',1,Nd,nSeries);
    switch SensMethod

        case 'brute-force'
            Nsim = N + Np2*N^2; % Number of Montecarlo simulations
            [Y] = gsua_pardeval(M,ParT,xdata,Parallel,[0 Nsim],verbose); % Y is a matrix (NxNd) with time responses in rows
            Nd = size(y_exp,2); % Number of simulation data
            J = sum((Y-yexp3).^2,[2 3]); % Scalar model output (summed over time and series)
            % Generation of (1xNp) cell with all combinations of parameters for calculation of
            % fractional variances. Every matrix has dimension (N^2xNp)
            % Mi{i} = [p1(1) ... pi(1) ... pNp(1);
            %          ...        ...       ...
            %          p1(N) ... pi(1) ... pNp(N);
            %          ---------------------------
            %          ...        ...       ...
            %          ---------------------------
            %          p1(1) ... pi(N) ... pNp(1);
            %          ...        ...       ...
            %          p1(N) ... pi(N) ... pNp(N)]
            Mi = cell(1,Np);
            for k = 1:Np
                Mi{k} = repmat(M,N,1); % Creation of a large matrix consisting of tiling of N copies of M
                for i = 1:N
                    Mi{k}((i-1)*N+1:i*N,k) = M(i,k)*ones(N,1); % The respective column has a constant value
                end
            end
            Yi = cell(1,Np);
            Ji = cell(1,Np);
            for k = 1:Np
                if ~(fixed(k))
                % Simulation with every set of parameters (row of Mi{k}) and storage
                % in matrix Yi (a cell for every parameter)
                [Yi{k}] = gsua_pardeval(Mi{k},ParT,xdata,Parallel,[N+(k-1)*N^2,Nsim],verbose);
                Ji{k} = sum((Yi{k}-yexp3).^2,[2 3]);
                end
            end
            V_vec = var(Y,0,1); % Total variance in every sample time of every column (1xNdxnSeries)
            Vi_vec = zeros(Np,Nd,nSeries);
            VJ = var(J,0);
            VJi = zeros(Np,1);
            for k = 1:Np
                if ~(fixed(k))
                ymean = zeros(N,Nd,nSeries);
                Jmean = zeros(N,1);
                for i = 1:N
                    ymean(i,:,:) = mean(Yi{k}((i-1)*N+1:i*N,:,:),1);
                    Jmean(i) = mean(Ji{k}((i-1)*N+1:i*N,:));
                end
                Vi_vec(k,:,:) = var(ymean,0,1); % Time-dependent fractional variances per series
                % Vi = V[E(Y|p(i)=p*)] with p*=[pi1 ... piN]'
                VJi(k) = var(Jmean,0);
                end
            end
            Si_vec = zeros(Np,Nd,nSeries);
            Si = zeros(Np,1);
            for k = 1:Np
                if ~(fixed(k))
                Si_vec(k,:,:) = Vi_vec(k,:,:)./V_vec; % Time-dependent fractional sensitivity indices
                Si(k) = VJi(k)/VJ;
                end
            end
            VTi_vec = zeros(Np,Nd,nSeries); % Time-dependent total variances: VTi_vec=E[V(Y|p(~i)=p*)]
            VTJi = zeros(Np,1);
            for k = 1:Np
                if ~(fixed(k))
                vary = zeros(N,Nd,nSeries);
                varJ = zeros(N,1);
                for i = 1:N
                    vary(i,:,:) = var(Yi{k}(i:N:(N-1)*N+i,:,:),0,1);
                    varJ(i) = var(Ji{k}(i:N:(N-1)*N+i,:),0);
                end
                VTi_vec(k,:,:) = mean(vary,1);
                VTJi(k) = mean(varJ);
                end
            end
            STi_vec = zeros(Np,Nd,nSeries);
            STi = zeros(Np,1);
            for k = 1:Np
                if ~(fixed(k))
                STi_vec(k,:,:) = VTi_vec(k,:,:)./V_vec; % Time-dependent total sensitivity
                STi(k) = VTJi(k)/VJ;
                end
            end

        case {'Sobol','Jansen','Saltelli','Xiao'}
            % A,B: Matrices (N/2xNp) with N/2 combinations of Np parameters. One parameter by column and
            % one sample by row. For example:
            % A = [p1(1)   p2(1)   ... pNp(1); p1(2)   p2(2) ... pNp(2); ... ; p1(N/2) p2(N/2) ... pNp(N/2)]
            A = M(1:N/2,:);
            B = M(N/2+1:N,:);
            % BAi matrices formed by all rows of B except the k column, which comes from A
            BAi = cell(1,Np);
            for k = 1:Np
                BAi{k} = B;
                BAi{k}(:,k) = A(:,k);
            end
            % ABi matrices formed by all colums of A except the k column, which comes from B
            ABi = cell(1,Np);
            for k = 1:Np
                ABi{k} = A;
                ABi{k}(:,k) = B(:,k);
            end

            if strcmp(SensMethod, 'Sobol')
                Nsim = N*(Np2+1); % Number of Montecarlo simulations
            else
                Nsim = N*(Np2/2+1); % Number of Montecarlo simulations
            end


            [YA] = gsua_pardeval(A,ParT,xdata,Parallel,[0,Nsim],verbose);
            [YB] = gsua_pardeval(B,ParT,xdata,Parallel,[N/2,Nsim],verbose); % dim(YA)=dim(YB)=(Nd x N/2)
            Y = [YA;YB];
            JA = sum((YA-yexp3).^2,[2 3]);
            JB = sum((YB-yexp3).^2,[2 3]);
            J = [JA;JB];
            Nd = size(y_exp,2); % Number of simulation data
            Si_vec = zeros(Np,Nd,nSeries);
            STi_vec = zeros(Np,Nd,nSeries);
            Si = zeros(Np,1);
            STi = zeros(Np,1);
            V_vec = var(Y,0,1); % Total variance (1xNd) of every column in every sample time
            VJ = var(J,0);
            YABi = cell(1,Np);
            YBAi = cell(1,Np);
            JABi = cell(1,Np);
            JBAi = cell(1,Np);

            switch SensMethod
                case 'Sobol'
                    f02_vec = mean(Y,1).^2;
                    f02 = mean(J)^2;
                    for k = 1:Np
                        if ~(fixed(k))
                            [YBAi{k}] = gsua_pardeval(BAi{k},ParT,xdata,Parallel,[k*N,Nsim],verbose);
                            [YABi{k}] = gsua_pardeval(ABi{k},ParT,xdata,Parallel,[k*N+N/2,Nsim],verbose);
                            JBAi{k} = sum((YBAi{k}-yexp3).^2,[2 3]);
                            JABi{k} = sum((YABi{k}-yexp3).^2,[2 3]);
                            Si_vec(k,:,:) = ( mean( YA.*YBAi{k},1 ) - f02_vec )./V_vec;
                            STi_vec(k,:,:) = mean( YA.*(YA - YABi{k}),1 )./V_vec;
                            Si(k) = ( mean( JA.*JBAi{k} ) - f02 )/VJ;
                            STi(k) = mean( JA.*(JA - JABi{k}) )/VJ;
                        end
                    end
                case 'Jansen'
                    for k = 1:Np
                        if ~(fixed(k))
                            [YABi{k}] = gsua_pardeval(ABi{k},ParT,xdata,Parallel,[N+(k-1)*N/2,Nsim],verbose);
                            JABi{k} = sum((YABi{k}-yexp3).^2,[2 3]);
                            Si_vec(k,:,:) = 1 - mean((YB - YABi{k}).^2,1)./(2*V_vec);
                            STi_vec(k,:,:) = mean( (YA - YABi{k}).^2,1 )./(2*V_vec);
                            Si(k) = 1 - mean((JB - JABi{k}).^2)/(2*VJ);
                            STi(k) = mean((JA - JABi{k}).^2)/(2*VJ);
                        end
                    end
                case 'Saltelli'
                    Y2=(Y-yexp3).^2;
                    V_vec2 = var(Y2,0,1);
                    YA2=(YA-yexp3).^2;
                    YB2=(YB-yexp3).^2;
                    for k = 1:Np
                        if ~(fixed(k))
                            [YABi{k}] = gsua_pardeval(ABi{k},ParT,xdata,Parallel,[N+(k-1)*N/2,Nsim],verbose);
                            JABi{k} = sum((YABi{k}-yexp3).^2,[2 3]);
                            YABi2=(YABi{k}-yexp3).^2;
                            Si_vec(k,:,:) = mean(YB2.*(YABi2 - YA2),1)./V_vec2;
                            STi_vec(k,:,:) = mean((YA2 - YABi2).^2,1)./(2*V_vec2);
                            %Si(k) = sum(mean(YB.*(YABi{k} - YA),1))/sum(V_vec);
                            %STi(k) = sum(mean((YA - YABi{k}).^2,1))/sum((2*V_vec));
                            Si(k) = mean(JB.*(JABi{k} - JA))/VJ;
                            STi(k) = mean((JA - JABi{k}).^2)/(2*VJ);
                        end
                    end
                case 'Xiao'
                    den=mean(vecnorm(reshape(YA-YB,size(YA,1),[]),2).^pod);  % scalar: norm over time & series
                    denv=mean((abs(YA-YB)).^pod);                            % per series: (1 x Nd x nSeries)
                        for k = 1:Np
                            if ~(fixed(k))
                                [YABi{k}] = gsua_pardeval(ABi{k},ParT,xdata,Parallel,[N+(k-1)*N/2,Nsim],verbose);
                                JABi{k} = sum((YABi{k}-yexp3).^2,[2 3]);
                                Si_vec(k,:,:)=(denv-mean((abs(YB-YABi{k})).^pod))./denv;
                                STi_vec(k,:,:) = mean((abs(YA-YABi{k})).^pod)./denv;
                                Si(k) = (den-mean(vecnorm(reshape(YB-YABi{k},size(YB,1),[]),2).^pod))/den;
                                STi(k) = mean(vecnorm(reshape(YA-YABi{k},size(YA,1),[]),2).^pod)/den;
                            end 
                        end
            end

        case 'Shapley'
            if isempty(shapPerm)
                shapPerm = N;
            end
            switch shapType
                case 'Independent'
                    [Sh,Sh_vec,Y,J] = local_shapley_indep(M,ParT,PT,xdata,yexp3,fixed, ...
                        Parallel,verbose,shapPerm,shapNo,shapNi);
                case 'Given'
                    if isempty(shapPool)
                        Xg = M;
                    elseif size(shapPool,1)==Np
                        Xg = shapPool.';        % Np x nPool ensemble (as Tia.Est) -> rows
                    elseif size(shapPool,2)==Np
                        Xg = shapPool;
                    else
                        error('gsua_sa:ShapleyPoolSize', ...
                            '''Pool'' must be Np x nPool (Np = %d parameters); got %d x %d.', ...
                            Np, size(shapPool,1), size(shapPool,2));
                    end
                    [Sh,Sh_vec,Y,J] = local_shapley_given(Xg,ParT,xdata,yexp3,fixed, ...
                        Parallel,verbose,shapK,shapPerm);
            end

        case 'OAT'
            % One factor at a time: factor i takes its N sampled values from M while every other
            % factor stays at T.Nominal -- the definition the Python port uses. Before v1.19 the
            % model received the single column M(:,i), so every evaluation failed.
            Nsim=N*Np2+N;
            [Y]=gsua_pardeval(M,ParT,xdata,Parallel,[0,Nsim],verbose);
            J = sum((Y-yexp3).^2,[2 3]);            % cost of the base sample, as for every method
            nominalRows = repmat(reshape(ParT.Nominal,1,[]),N,1);
            Voat = zeros(Np,1);
            done = 0;
            for i = 1:Np
                if ~(fixed(i))
                    M_oat = nominalRows;
                    M_oat(:,i) = M(:,i);
                    done = done + 1;
                    [Y_oat]=gsua_pardeval(M_oat,ParT,xdata,Parallel,[done*N,Nsim],verbose);
                    Voat(i) = var(sum((Y_oat-yexp3).^2,[2 3]),1);   % population variance (= Python)
                end
            end
            VT = sum(Voat);
            if VT > 0
                Si = Voat/VT;
            else
                Si = zeros(Np,1);
            end
            STi = Voat;

        otherwise
            clc
            disp('---------------------------------------')
            disp('Error ... Give the right name of Sensitivity Method')
            disp('---------------------------------------')
            return
    end
    if strcmp(SensMethod,'Shapley')
        ParT.Sh=Sh;
        ParT.Sh_vec=Sh_vec;
    elseif strcmp(SensMethod,'OAT')
        ParT.Si=Si;                 % OAT has no time-dependent indices (as in Python)
        ParT.STi=STi;
    else
        ParT.Si=Si;
        ParT.STi=STi;
        ParT.STi_vec=STi_vec;
        ParT.Si_vec=Si_vec;
    end

else %Multiobjetive SA
    if strcmp(SensMethod,'Shapley')
        error('gsua_sa:ShapleyMultiObjective', ['''Shapley'' is not available with ' ...
            '''multiObjective'',true; run it per output series (the default).']);
    end
    switch SensMethod
    
        case 'brute-force'
            Nsim = N + Np2*N^2; % Number of Montecarlo simulations
            [Y] = gsua_pardeval(M,ParT,xdata,Parallel,[0 Nsim],verbose); % Y is a matrix (NxNd) with time responses in rows
            if isempty(y_exp)
                disp('Finding the deepest curve...')
                [~,idx]=gsua_depth(Y,Parallel);
                y_exp=Y(idx(1),:,:);
            end
            Y = gsua_costfMulti(Y,y_exp,margin,Parallel); % Scalar model output
            % Generation of (1xNp) cell with all combinations of parameters for calculation of
            % fractional variances. Every matrix has dimension (N^2xNp)
            % Mi{i} = [p1(1) ... pi(1) ... pNp(1);
            %          ...        ...       ...
            %          p1(N) ... pi(1) ... pNp(N);
            %          ---------------------------
            %          ...        ...       ...
            %          ---------------------------
            %          p1(1) ... pi(N) ... pNp(1);
            %          ...        ...       ...
            %          p1(N) ... pi(N) ... pNp(N)]
            Mi = cell(1,Np);
            for k = 1:Np
                Mi{k} = repmat(M,N,1); % Creation of a large matrix consisting of tiling of N copies of M
                for i = 1:N
                    Mi{k}((i-1)*N+1:i*N,k) = M(i,k)*ones(N,1); % The respective column has a constant value
                end
            end
            Yi = cell(1,Np);
            for k = 1:Np
                if ~(fixed(k))
                % Simulation with every set of parameters (row of Mi{k}) and storage
                % in matrix Yi (a cell for every parameter)
                [Yi{k}] = gsua_costfMulti(gsua_pardeval(Mi{k},ParT,xdata,Parallel,[N+(k-1)*N^2,Nsim],verbose),...
                    y_exp,margin,Parallel);
                end
            end
            VJ = var(Y,0);
            VJi = zeros(Np,1);
            for k = 1:Np
                if ~(fixed(k))
                ymean = zeros(N,1);
                for i = 1:N
                    ymean(i,:) = mean(Yi{k}((i-1)*N+1:i*N,:),1);
                end
                VJi(k) = var(ymean,0);
                end
            end
            Si = zeros(Np,1);
            for k = 1:Np
                if ~(fixed(k))
                    Si(k) = VJi(k)/VJ;
                end
            end
            VTJi = zeros(Np,1);
            for k = 1:Np
                if ~(fixed(k))
                vary = zeros(N,1);
                for i = 1:N
                    vary(i,:) = var(Yi{k}(i:N:(N-1)*N+i,:),0,1);
                end
                VTJi(k) = mean(vary);
                end
            end
            STi = zeros(Np,1);
            for k = 1:Np
                if ~(fixed(k))
                    STi(k) = VTJi(k)/VJ;
                end
            end

        case {'Sobol','Jansen','Saltelli','Xiao'}
            % A,B: Matrices (N/2xNp) with N/2 combinations of Np parameters. One parameter by column and
            % one sample by row. For example:
            % A = [p1(1)   p2(1)   ... pNp(1); p1(2)   p2(2) ... pNp(2); ... ; p1(N/2) p2(N/2) ... pNp(N/2)]
            A = M(1:N/2,:);
            B = M(N/2+1:N,:);
            % BAi matrices formed by all rows of B except the k column, which comes from A
            BAi = cell(1,Np);
            for k = 1:Np
                BAi{k} = B;
                BAi{k}(:,k) = A(:,k);
            end
            % ABi matrices formed by all colums of A except the k column, which comes from B
            ABi = cell(1,Np);
            for k = 1:Np
                ABi{k} = A;
                ABi{k}(:,k) = B(:,k);
            end

            if strcmp(SensMethod, 'Sobol')
                Nsim = N*(Np2+1); % Number of Montecarlo simulations
            else
                Nsim = N*(Np2/2+1); % Number of Montecarlo simulations
            end


            [YA] = gsua_pardeval(A,ParT,xdata,Parallel,[0,Nsim],verbose);
            [YB] = gsua_pardeval(B,ParT,xdata,Parallel,[N/2,Nsim],verbose); % dim(YA)=dim(YB)=(Nd x N/2)
            Y = [YA;YB];
            if isempty(y_exp)
                disp('Finding the deepest curve...')
                [~,idx]=gsua_depth(Y,Parallel);
                y_exp=Y(idx(1),:,:);
            end
            YA = gsua_costfMulti(YA,y_exp,margin,Parallel); % Scalar model output
            YB = gsua_costfMulti(YB,y_exp,margin,Parallel); % Scalar model output
            Y = [YA,YB];
            Si = zeros(Np,1);
            STi = zeros(Np,1);
            VJ = var(Y,0);
            YABi = cell(1,Np);
            YBAi = cell(1,Np);

            switch SensMethod
                case 'Sobol'
                    f02 = mean(Y)^2;
                    for k = 1:Np
                        if ~(fixed(k))
                            [YBAi{k}] = gsua_costfMulti(gsua_pardeval(BAi{k},ParT,xdata,Parallel,[k*N,Nsim],verbose),...
                                y_exp,margin,Parallel);
                            [YABi{k}] = gsua_costfMulti(gsua_pardeval(ABi{k},ParT,xdata,Parallel,[k*N+N/2,Nsim],verbose),...
                                y_exp,margin,Parallel);
                            Si(k) = ( mean( YA.*YBAi{k} ) - f02 )/VJ;
                            STi(k) = mean( YA.*(YA - YABi{k}) )/VJ;
                        end
                    end
                case 'Jansen'
                    for k = 1:Np
                        if ~(fixed(k))
                            [YABi{k}] = gsua_costfMulti(gsua_pardeval(ABi{k},ParT,xdata,Parallel,[N+(k-1)*N/2,Nsim],verbose),...
                                y_exp,margin,Parallel);
                            Si(k) = 1 - mean((YB - YABi{k}).^2)/(2*VJ);
                            STi(k) = mean((YA - YABi{k}).^2)/(2*VJ);
                        end
                    end
                case 'Saltelli'
                    for k = 1:Np
                        if ~(fixed(k))
                            [YABi{k}] = gsua_costfMulti(gsua_pardeval(ABi{k},ParT,xdata,Parallel,[N+(k-1)*N/2,Nsim],verbose),...
                                y_exp,margin,Parallel);
                            Si(k) = mean(YB.*(YABi{k} - YA))/VJ;
                            STi(k) = mean((YA - YABi{k}).^2)/(2*VJ);
                        end
                    end
                case 'Xiao'
                    den=mean(vecnorm(YA-YB,2).^pod);
                        for k = 1:Np
                            if ~(fixed(k))
                                [YABi{k}] = gsua_costfMulti(gsua_pardeval(ABi{k},ParT,xdata,Parallel,[N+(k-1)*N/2,Nsim],verbose),...
                                    y_exp,margin,Parallel);
                                Si(k) = (den-mean(vecnorm(YB-YABi{k},2).^pod))/den;
                                STi(k) = mean(vecnorm(YA-YABi{k},2).^pod)/den;
                            end 
                        end
            end

            otherwise
                clc
                disp('---------------------------------------')
                disp('Error ... Not supported SA method: OAT')
                disp('---------------------------------------')
            return
    end
    ParT.Si=Si;
    ParT.STi=STi;

end
try
    try
    ParT=addprop(ParT,{'SensMethod'},{'table'});
    catch ME
        warning(ME.message)
    end
    ParT.Properties.CustomProperties.SensMethod=SensMethod;
catch
    Table=load('ATable.mat');
    Table2=Table.Table2;
    Table2.SensMethod=SensMethod;
    save('ATable','Table2');
end

end

% =====================================================================================================
% Shapley effects
% Value of a coalition J (the factors that are KNOWN): c(J) = E[Var(Y | X_J)], so c({}) = Var(Y) and
% c(all) = 0 for a deterministic model. Sh_i = (1/Var Y) * E_perm[ c(P_i) - c(P_i U {i}) ], with P_i the
% factors before i in a uniformly random permutation. The increments telescope from Var(Y) to 0 along
% every permutation, so the indices sum to 1 exactly. The scalar index uses the cost
% J = sum((Y-ynom).^2) as the output (like Si/STi); the vector index uses Y at every time and series.
% =====================================================================================================

function [Sh,Sh_vec,Y,J] = local_shapley_indep(M,ParT,PT,xdata,yexp3,fixed,Parallel,verbose,m,No,Ni)
% Song, Nelson & Staum (2016) nested Monte Carlo for independent factors: for each prefix coalition
% of a permutation, draw No outer values of the known factors and Ni inner values of the others from
% their marginals, and average the inner variances.
[N,Np] = size(M);
Nd = size(yexp3,2);
nSeries = size(yexp3,3);
fixed = logical(fixed(:).');
free = find(~fixed);
d = numel(free);
Nsim = N + m*max(d-1,0)*No*Ni;

Y = gsua_pardeval(M,ParT,xdata,Parallel,[0,Nsim],verbose);
J = sum((Y-yexp3).^2,[2 3]);
Vv = reshape(var(Y,0,1),1,Nd*nSeries);
VJ = var(J,0);
Sh = zeros(Np,1);
Sh_vec = zeros(Np,Nd*nSeries);
if d==1
    Sh(free) = 1;
    Sh_vec(free,:) = Vv./Vv;
end
if d<=1
    Sh_vec = reshape(Sh_vec,Np,Nd,nSeries);
    return
end

% Permutations of the free factors. When every permutation fits in the budget, use each one equally
% often (stratified) -- lower variance than drawing them at random.
nP = factorial(d);
if nP <= m
    Pall = perms(1:d);
    reps = floor(m/nP);
    P = [repmat(Pall,reps,1); Pall(randperm(nP,m-reps*nP),:)];
else
    P = zeros(m,d);
    for i = 1:m
        P(i,:) = randperm(d);
    end
end

% Group g = ((i-1)*(d-1) + j-1)*No + o holds Ni rows that share the outer draw of the coalition
% free(P(i,1:j)) and differ in every other free factor.
G = m*(d-1)*No;
known = false(G,Np);
g = 0;
for i = 1:m
    for j = 1:d-1
        row = false(1,Np);
        row(free(P(i,1:j))) = true;
        known(g+1:g+No,:) = repmat(row,No,1);
        g = g + No;
    end
end
outer = local_marginal_sample(M,ParT,PT,fixed,G);
X = local_marginal_sample(M,ParT,PT,fixed,G*Ni);
grp = repelem((1:G).',Ni);
K = known(grp,:);
O = outer(grp,:);
X(K) = O(K);

Yb = gsua_pardeval(X,ParT,xdata,Parallel,[N,Nsim],verbose);
Jb = sum((Yb-yexp3).^2,[2 3]);
cv = var(reshape(Yb,Ni,G,Nd*nSeries),0,1);                          % 1 x G x (Nd*nSeries)
cv = reshape(mean(reshape(cv,No,m*(d-1),Nd*nSeries),1),d-1,m,Nd*nSeries);
cJ = reshape(mean(reshape(var(reshape(Jb,Ni,G),0,1),No,m*(d-1)),1),d-1,m);

acc = zeros(Np,Nd*nSeries);
accJ = zeros(Np,1);
zeroRow = zeros(1,Nd*nSeries);
for i = 1:m
    prev = Vv;
    prevJ = VJ;
    for j = 1:d
        if j<d
            C = reshape(cv(j,i,:),1,[]);
            CJ = cJ(j,i);
        else
            C = zeroRow;
            CJ = 0;
        end
        k = free(P(i,j));
        acc(k,:) = acc(k,:) + (prev - C);
        accJ(k) = accJ(k) + (prevJ - CJ);
        prev = C;
        prevJ = CJ;
    end
end
Sh = accJ/(m*VJ);
Sh_vec = acc./(m*Vv);
Sh(fixed) = 0;
Sh_vec(fixed,:) = 0;
Sh_vec = reshape(Sh_vec,Np,Nd,nSeries);
end

function [Sh,Sh_vec,Y,J] = local_shapley_given(X,ParT,xdata,yexp3,fixed,Parallel,verbose,k,m)
% Broto, Bachoc & Depecker (2020) given-data estimator: one joint sample, c(J) estimated by the
% variance of Y over the k nearest neighbours of each point in the (standardized) X_J coordinates.
% Exact Shapley formula over all coalitions up to 12 players, random permutations above.
X = unique(X,'rows','stable');      % a repeated vector adds no information and would zero the kNN variance
[n,Np] = size(X);
fixed = logical(fixed(:).');
Nd = size(yexp3,2);
nSeries = size(yexp3,3);
if n <= k
    error('gsua_sa:ShapleySampleSize', ...
        'The joint sample has %d distinct vectors; ''Given'' needs more than ''ShapleyK''=%d.', n, k);
end
if n < 100
    warning('gsua_sa:ShapleySmallSample', ['Only %d distinct vectors in the joint sample: the ' ...
        'nearest-neighbour Shapley estimates will be noisy. Use a few hundred or more.'], n);
end
Y = gsua_pardeval(X,ParT,xdata,Parallel,[0,n],verbose);
J = sum((Y-yexp3).^2,[2 3]);
Yf = reshape(Y,n,Nd*nSeries);
Vv = var(Yf,0,1);
VJ = var(J,0);
sd = std(X,0,1);
players = find(~fixed(:).' & sd>0);       % a factor constant in the sample cannot carry variance
d = numel(players);
Sh = zeros(Np,1);
Sh_vec = zeros(Np,Nd*nSeries);
if d==1
    Sh(players) = 1;
    Sh_vec(players,:) = Vv./Vv;
end
if d<=1
    Sh_vec = reshape(Sh_vec,Np,Nd,nSeries);
    return
end
Z = (X(:,players) - mean(X(:,players),1))./sd(players);
zeroRow = zeros(1,Nd*nSeries);

if d <= 12
    nS = 2^d;
    C = zeros(nS,Nd*nSeries);
    CJ = zeros(nS,1);
    C(1,:) = Vv;                      % empty coalition; the full one (row nS) stays 0
    CJ(1) = VJ;
    for mask = 1:nS-2
        [C(mask+1,:),CJ(mask+1)] = local_knnvar(Z(:,bitget(mask,1:d)==1),Yf,J,k);
    end
    for a = 1:d
        for mask = 0:nS-1
            if bitget(mask,a)
                continue
            end
            w = 1/(d*nchoosek(d-1,sum(bitget(mask,1:d))));
            mi = bitset(mask,a);
            Sh_vec(players(a),:) = Sh_vec(players(a),:) + w*(C(mask+1,:) - C(mi+1,:));
            Sh(players(a)) = Sh(players(a)) + w*(CJ(mask+1) - CJ(mi+1));
        end
    end
else
    memo = containers.Map('KeyType','double','ValueType','any');
    for i = 1:m
        pp = randperm(d);
        mask = 0;
        prev = Vv;
        prevJ = VJ;
        for j = 1:d
            mask = bitset(mask,pp(j));
            if j==d
                Cc = zeroRow;
                CJc = 0;
            elseif isKey(memo,mask)
                v = memo(mask);
                Cc = v{1};
                CJc = v{2};
            else
                [Cc,CJc] = local_knnvar(Z(:,bitget(mask,1:d)==1),Yf,J,k);
                memo(mask) = {Cc,CJc};
            end
            Sh_vec(players(pp(j)),:) = Sh_vec(players(pp(j)),:) + (prev - Cc);
            Sh(players(pp(j))) = Sh(players(pp(j))) + (prevJ - CJc);
            prev = Cc;
            prevJ = CJc;
        end
    end
    Sh_vec = Sh_vec/m;
    Sh = Sh/m;
end
Sh = Sh/VJ;
Sh_vec = Sh_vec./Vv;
nonPlayers = true(1,Np);
nonPlayers(players) = false;
Sh(nonPlayers) = 0;
Sh_vec(nonPlayers,:) = 0;
Sh_vec = reshape(Sh_vec,Np,Nd,nSeries);
end

function [C,CJ] = local_knnvar(Zj,Yf,J,k)
% Mean over the sample of the variance of the output among each point's k nearest neighbours in Zj.
n = size(Zj,1);
idx = knnsearch(Zj,Zj,'K',k);
C = reshape(mean(var(reshape(Yf(idx(:),:),n,k,[]),0,2),1),1,[]);
CJ = mean(var(J(idx),0,2));
end

function X = local_marginal_sample(M,ParT,PT,fixed,n)
% Independent draws from each factor's marginal, matching gsua_dmatrix: N(mu,sigma) when the ranges
% are normal parameters (rMethod 'normal'), uniform on T.Range otherwise. Fixed factors keep the
% value they have in M.
Np = size(M,2);
R = ParT.Range;
isNormal = false;
try
    isNormal = any(strcmp(PT.rMethod,{'normal'}));
catch
end
X = zeros(n,Np);
for k = 1:Np
    if fixed(k)
        X(:,k) = M(1,k);
    elseif isNormal
        X(:,k) = R(k,1) + R(k,2)*randn(n,1);
    else
        X(:,k) = R(k,1) + (R(k,2)-R(k,1))*rand(n,1);
    end
end
end