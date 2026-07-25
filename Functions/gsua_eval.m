function [y,xdata] = gsua_eval(par,Table,varargin)
% Function for few model evaluations
%
% Usage:
%   Y = gsua_eval(values,Table)
%   Y = gsua_eval(values,Table,xdata)
%   Y = gsua_eval(values,Table,xdata,ydata)
%   Y = gsua_eval(values,Table,xdata,ydata,parallel)
%   Y = gsua_eval(values,Table,xdata,ydata,parallel,show)
%
% Inputs:
%   par      <-- array of Np x N (number of factors x number of simulations)
%   Table    <-- summary table from gsua_dataprep function
%   xdata    <-- (optional) domain to which model output is interpolated
%   ydata    <-- (optional) experimental / previous model data matching xdata
%   parallel <-- (optional, logical) use parallel evaluation (default: false)
%   show     <-- (optional, logical) if true, calls gsua_plot; if false, no plot
%
% Outputs:
%   y     <-- array with model output
%   xdata <-- domain used for evaluation/interpolation
%
% See also GSUA_DEVAL, GSUA_PARDEVAL, GSUA_PLOT.

    % --------- Parse inputs smoothly ---------
    ip = inputParser;
    ip.FunctionName = 'gsua_eval';

    addRequired(ip,'par');
    addRequired(ip,'Table');
    addOptional(ip,'xdata',   [], @(x) true);
    addOptional(ip,'ydata',   [], @(x) true);
    addOptional(ip,'parallel',false,@(x) islogical(x) || isnumeric(x));
    addOptional(ip,'show',    true, @(x) islogical(x) || isnumeric(x));

    parse(ip,par,Table,varargin{:});

    par      = ip.Results.par;
    Table    = ip.Results.Table;
    xdata    = ip.Results.xdata;
    ydata    = ip.Results.ydata;
    parallel = logical(ip.Results.parallel);
    show     = logical(ip.Results.show);

    % --------- Parallel configuration ---------
    if parallel
        parforArg = inf;
    else
        parforArg = 0;
    end

    % Transpose parameter matrix as in original code
    par = par';

    % --------- Recover custom properties from Table ---------
    try
        Table2 = Table.Properties.CustomProperties;
    catch
        Table2 = load('ATable.mat');
        Table2 = Table2.Table2;
    end

    kind = Table2.Kind;

    % --------- Default xdata if not provided ---------
    if isempty(xdata)
        try
            xdata = Table2.Domain(1):Table2.Domain(2);
        catch
            xdata = 1;
        end
    end

    % --------- Main evaluation logic ---------
    if kind == 1
        % Possibly it is wrong (kept as in original code)
        [y,~] = sens_montecarlo(Table,par,xdata,false,[0,size(par,1)]);
    else
        fun  = Table2.Solver;
        ni   = Table2.NumVars;
        out  = Table2.output;
        nout = size(out,2);

        % total number of simulations
        nSim = size(par,1);
        % number of blocks (at most 10, but not more than nSim)
        nBlocks = min(10, nSim);
        % block edges (indices)
        edges = round(linspace(1, nSim+1, nBlocks+1));

        switch kind
            case 3
                y = zeros(nSim, numel(xdata), nout);

                % special case: exactly one simulation (keep old behavior with y_temp)
                if nSim == 1
                    y_temp = deval(fun(par(1,1:ni),par(1,ni+1:end)),xdata,out);
                    for j = 1:nout
                        y(1,:,j) = y_temp(j,:);
                    end
                    % progress is trivially 100%
                    disp('Progress: 100 %');
                else
                    % multi-simulation case in blocks
                    for b = 1:nBlocks
                        iStart = edges(b);
                        iEnd   = edges(b+1) - 1;

                        if parallel
                            parfor (i = iStart:iEnd, parforArg)
                                y_temp2 = deval(fun(par(i,1:ni),par(i,ni+1:end)),xdata,out);
                                for j = 1:nout
                                    y(i,:,j) = y_temp2(j,:);
                                end
                            end
                        else
                            for i = iStart:iEnd
                                y_temp2 = deval(fun(par(i,1:ni),par(i,ni+1:end)),xdata,out);
                                for j = 1:nout
                                    y(i,:,j) = y_temp2(j,:);
                                end
                            end
                        end

                        progress = 100 * iEnd / nSim;
                        disp(['Progress: ', num2str(progress,'%.1f'), ' %']);
                    end
                end

            case 2
                % black-box vectorized solver, no internal loop here
                y = fun(par);

            case 4
                % black-box vectorized solver as well
                y = gsua_intrp(fun(par(:,1:ni),par(:,ni+1:end)),xdata,out);

            case 5
                domain = Table2.Domain;
                opt    = Table2.copt;
                y      = zeros(nSim, numel(xdata), nout);

                for b = 1:nBlocks
                    iStart = edges(b);
                    iEnd   = edges(b+1) - 1;

                    if parallel
                        parfor (i = iStart:iEnd, parforArg)
                            y(i,:,:) = gsua_intrp(fun(par(i,:),domain,opt),xdata,out);
                        end
                    else
                        for i = iStart:iEnd
                            y(i,:,:) = gsua_intrp(fun(par(i,:),domain,opt),xdata,out);
                        end
                    end

                    progress = 100 * iEnd / nSim;
                    disp(['Progress: ', num2str(progress,'%.1f'), ' %']);
                end

            case 6
                y = zeros(nSim,1);

                for b = 1:nBlocks
                    iStart = edges(b);
                    iEnd   = edges(b+1) - 1;

                    if parallel
                        parfor (i = iStart:iEnd, parforArg)
                            y(i) = fun(par(i,:));
                        end
                    else
                        for i = iStart:iEnd
                            y(i) = fun(par(i,:));
                        end
                    end

                    progress = 100 * iEnd / nSim;
                    disp(['Progress: ', num2str(progress,'%.1f'), ' %']);
                end
        end
    end

    % --------- Optional plotting (controlled by "show") ---------
    if show
        if ~isempty(ydata)
            gsua_plot('UncertaintyAnalysis',Table,y,xdata,ydata);
        else
            gsua_plot('UncertaintyAnalysis',Table,y,xdata);
        end
    end

    % --------- Keep original special case behavior ---------
    if size(y,1) == 1 && (kind == 3)
        % y_temp exists only in the nSim == 1 branch above
        y = y_temp;
    end

end
