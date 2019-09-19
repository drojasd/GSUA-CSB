function [ParT,J,Y] = gsua_sa(M,ParT,varargin)
% Function for sensitivity analysis
%
% T=gsua_sa(M,T)
% Parameters:
% M <-- design matrix from gsua_dmatrix function
% T <-- summary table from gsua_dataprep function
% Output:
% T <-- summary table with sensitivity indices.
% Additional features:
% It is possible to request 2 extra positional outputs:
% J <-- vector of model outputs in scalar representation (MSE)
% Y <-- model output matrix (same as gsua_ua output)
% It is possible to choose between 6 SA methods using the paired feature
% 'SensMethod',method. Where method is one of the following:
% 'Saltelli'    <-- requires N*(Np/2+1) simulations
% 'Jansen'      <-- requires N*(Np/2+1) simulations
% 'Xiao'        <-- newest method (Default), requires N*(Np/2+1)
%                   simulations
% 'Sobol'       <-- classic method, requires N*(Np+1) simulations
% 'brute-force' <-- explores all possible combinations, N + Np2*N^2
%                   simulations
% 'OAT'         <-- non-global method, requires N*(Np+1) simulations
% Also, it is possible to specify an output for indices estimation using
% the paired feature 'ynom'. Finally, you can avoid parallel speed up 
% for simulation process using the paired feature 'parallel', false.
% [ParT,J,Y] = gsua_sa(M,T,'SensMethod','Xiao','ynom',ynom)
% Note: ynom size must coincide with xdata size, where xdata is the next
% xdata=linspace(Domain(1),Domain(2),Domain(2)-Domain(1)+1) and Domain is
% the array of the model domain


p=inputParser;

defaultSensMethod='Xiao';
validSensMethod={'brute-force','Jansen','Sobol','OAT','Saltelli','Xiao'};
checkSensMethod = @(x) any(validatestring(x,validSensMethod));

defaultParallel=true;

defaultSim_progress='off';
checkSim_progress = @(x) any(validatestring(x,validParallel));

defaultyexp=[];


addRequired(p,'M',@isnumeric);
addRequired(p,'ParT',@istable);
addParameter(p,'ynom',defaultyexp,@isnumeric);
addParameter(p,'SensMethod',defaultSensMethod,checkSensMethod);
addParameter(p,'parallel',defaultParallel,@islogical);
addParameter(p,'Sim_progress',defaultSim_progress,checkSim_progress);

parse(p,M,ParT,varargin{:})

M=p.Results.M;
ParT=p.Results.ParT;
y_exp=p.Results.ynom;
SensMethod=p.Results.SensMethod;
Parallel=p.Results.parallel;
try
    PT=ParT.Properties.CustomProperties;
catch 
    PT=load('ATable.mat');
    PT=PT.Table2;
end

Sim_progress=p.Results.Sim_progress;
kind=PT.Kind;
if strcmp('mat',kind)
domain=PT.Domain;
Sname=strcmp(PT.Sname,'ode45');
if size(domain,2)==1
    xdata=domain;
    Parallel=false;
else
    xdata=linspace(domain(1),domain(2),domain(2)-domain(1)+1);
end
else
    xdata=[];
    Sname=true;
end
Par=ParT.Properties.RowNames;
fixed=PT.Fixed;
if isempty(y_exp)
    y_exp=gsua_deval(ParT.Nominal',ParT,xdata);
end


if Parallel && exist('parpool','file')>0 && isempty(gcp('nocreate')) && Sname
    parpool local;
elseif ~Parallel && exist('parpool','file')>0 && isempty(gcp('nocreate'))==0
    delete(gcp('nocreate'))
end

Np=size(Par,1);
Np2=size(Par,1)-sum(fixed);
N=size(M,1);
if isempty(fixed)
    fixed=false(1,Np);
end

switch SensMethod
    
    case 'brute-force'
        Nsim = N + Np2*N^2; % Number of Montecarlo simulations
        [Y] = gsua_pardeval(M,ParT,xdata,Parallel,[0 Nsim]); % Y is a matrix (NxNd) with time responses in rows
        Nd = size(y_exp,2); % Number of simulation data
        Y_nom = ones(N,1)*y_exp;
        J = sum((Y-Y_nom).^2,2); % Scalar model output
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
            [Yi{k}] = gsua_pardeval(Mi{k},ParT,xdata,Parallel,[N+(k-1)*N^2,Nsim]);
            Y_nom2 = ones(N^2,1)*y_exp;
            Ji{k} = sum((Yi{k}-Y_nom2).^2,2);
            end
        end
        V_vec = var(Y,0,1); % Total variance in every sample time of every column
        Vi_vec = zeros(Np,Nd);
        VJ = var(J,0);
        VJi = zeros(Np,1);
        for k = 1:Np
            if ~(fixed(k))
            ymean = zeros(N,Nd);
            Jmean = zeros(N,1);
            for i = 1:N
                ymean(i,:) = mean(Yi{k}((i-1)*N+1:i*N,:),1);
                Jmean(i) = mean(Ji{k}((i-1)*N+1:i*N,:));
            end
            Vi_vec(k,:) = var(ymean,0,1); % Time-dependent fractional variances of every column
            % Vi = V[E(Y|p(i)=p*)] with p*=[pi1 ... piN]'
            VJi(k) = var(Jmean,0);
            end
        end
        Si_vec = zeros(Np,Nd);
        Si = zeros(Np,1);
        for k = 1:Np
            if ~(fixed(k))
            Si_vec(k,:) = Vi_vec(k,:)./V_vec; % Time-dependent fractional sensitivity indices
            Si(k) = VJi(k)/VJ;
            end
        end
        VTi_vec = zeros(Np,Nd); % Time-dependent total variances: VTi_vec=E[V(Y|p(~i)=p*)]
        VTJi = zeros(Np,1);
        for k = 1:Np
            if ~(fixed(k))
            vary = zeros(N,Nd);
            varJ = zeros(N,1);
            for i = 1:N
                vary(i,:) = var(Yi{k}(i:N:(N-1)*N+i,:),0,1);
                varJ(i) = var(Ji{k}(i:N:(N-1)*N+i,:),0);
            end
            VTi_vec(k,:) = mean(vary,1);
            VTJi(k) = mean(varJ);
            end
        end
        STi_vec = zeros(Np,Nd);
        STi = zeros(Np,1);
        for k = 1:Np
            if ~(fixed(k))
            STi_vec(k,:) = VTi_vec(k,:)./V_vec; % Time-depndent total sensitivity
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
        
        
        [YA] = gsua_pardeval(A,ParT,xdata,Parallel,[0,Nsim]);
        [YB] = gsua_pardeval(B,ParT,xdata,Parallel,[N/2,Nsim]); % dim(YA)=dim(YB)=(Nd x N/2)
        Y = [YA;YB];
        JA = sum((YA-ones(N/2,1)*y_exp).^2,2);
        JB = sum((YB-ones(N/2,1)*y_exp).^2,2);
        J = [JA;JB];
        Nd = size(y_exp,2); % Number of simulation data
        Si_vec = zeros(Np,Nd);
        STi_vec = zeros(Np,Nd);
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
                    [YBAi{k}] = gsua_pardeval(BAi{k},ParT,xdata,Parallel,[k*N,Nsim]);
                    [YABi{k}] = gsua_pardeval(ABi{k},ParT,xdata,Parallel,[k*N+N/2,Nsim]);
                    JBAi{k} = sum((YBAi{k}-ones(N/2,1)*y_exp).^2,2);
                    JABi{k} = sum((YABi{k}-ones(N/2,1)*y_exp).^2,2);
                    Si_vec(k,:) = ( mean( YA.*YBAi{k},1 ) - f02_vec )./V_vec;
                    STi_vec(k,:) = mean( YA.*(YA - YABi{k}),1 )./V_vec;
                    Si(k) = ( mean( JA.*JBAi{k} ) - f02 )/VJ;
                    STi(k) = mean( JA.*(JA - JABi{k}) )/VJ;
                    end
                end
            case 'Jansen'
                for k = 1:Np
                    if ~(fixed(k))
                    [YABi{k}] = gsua_pardeval(ABi{k},ParT,xdata,Parallel,[N+(k-1)*N/2,Nsim]);
                    JABi{k} = sum((YABi{k}-ones(N/2,1)*y_exp).^2,2);
                    Si_vec(k,:) = 1 - mean((YB - YABi{k}).^2,1)./(2*V_vec);
                    STi_vec(k,:) = mean( (YA - YABi{k}).^2,1 )./(2*V_vec);
                    Si(k) = 1 - mean((JB - JABi{k}).^2)/(2*VJ);
                    STi(k) = mean((JA - JABi{k}).^2)/(2*VJ);
                    end
                end
            case 'Saltelli'
                Y2=(Y-y_exp).^2;
                V_vec2 = var(Y2,0,1);
                YA2=(YA-y_exp).^2;
                YB2=(YB-y_exp).^2;
                for k = 1:Np
                    if ~(fixed(k))
                    [YABi{k}] = gsua_pardeval(ABi{k},ParT,xdata,Parallel,[N+(k-1)*N/2,Nsim]);
                    JABi{k} = sum((YABi{k}-ones(N/2,1)*y_exp).^2,2);
                    YABi2=(YABi{k}-y_exp).^2;
                    Si_vec(k,:) = mean(YB2.*(YABi2 - YA2),1)./V_vec2;
                    STi_vec(k,:) = mean((YA2 - YABi2).^2,1)./(2*V_vec2);
                    %Si(k) = sum(mean(YB.*(YABi{k} - YA),1))/sum(V_vec);
                    %STi(k) = sum(mean((YA - YABi{k}).^2,1))/sum((2*V_vec));
                    Si(k) = mean(JB.*(JABi{k} - JA))/VJ;
                    STi(k) = mean((JA - JABi{k}).^2)/(2*VJ);
                    end
                end
            case 'Xiao'
                pod=2;
                den=mean(vecnorm(YA-YB,2).^pod);
                    for k = 1:Np
                    if ~(fixed(k))
                    [YABi{k}] = gsua_pardeval(ABi{k},ParT,xdata,Parallel,[N+(k-1)*N/2,Nsim]);
                    JABi{k} = sum((YABi{k}-ones(N/2,1)*y_exp).^2,2);
                    Si_vec(k,:) = mean(YB.*(YABi{k} - YA),1)./V_vec;
                    STi_vec(k,:) = mean((YA - YABi{k}).^2,1)./(2*V_vec);
                    Si(k) = (den-mean(vecnorm(YB-YABi{k},2).^pod))/den;
                    STi(k) = mean(vecnorm(YA-YABi{k},2).^pod)/den;
                    end 
                    end
        end
        
    case 'OAT'
        if nargout>4
            disp('Too many output arguments')
            return
        else
            Nsim=N*Np2+N;
            [Y]=gsua_pardeval(M,ParT,xdata,Parallel,[0,Nsim]);
            J = zeros(Np,1);
            for i = 1:Np %construir M apropiado
                if ~(fixed(k))
                M_oat=M(:,i);
                [Y_oat]=gsua_pardeval(M_oat,ParT,xdata,Parallel,[i*N,Nsim]);
                J(i) = var(sum((Y_oat-y_exp).^2,2));
                end
            end

            VT = sum(J);
            Si = J/VT;
            STi=J;
        end
        
    otherwise
        clc
        disp('---------------------------------------')
        disp('Error ... Give the right name of Sensitivity Method')
        disp('---------------------------------------')
        return
end
try
    ParT=addprop(ParT,{'SensMethod'},{'table'});
    ParT.Properties.CustomProperties.SensMethod=SensMethod;
catch
    Table=load('ATable.mat');
    Table2=Table.Table2;
    Table2.SensMethod=SensMethod;
    save('ATable','Table2');
end
ParT.Si=Si;
ParT.STi=STi;
ParT.STi_vec=STi_vec;
ParT.Si_vec=Si_vec;


end