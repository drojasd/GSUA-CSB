function y = gsua_deval(par,Table2,xdata)
%GSUA_DEVAL Evaluate the model behind a GSUA table for one parameter set
%
%   y = gsua_deval(par,T,xdata)
%
%   Internal dispatcher used throughout the toolbox (GSUA_PE, GSUA_OATR,
%   GSUA_LIKELIHOOD, ...) to evaluate a single candidate parameter
%   vector against whatever model kind T was built from (Simulink,
%   symbolic ODE, or user-defined function) and return output
%   interpolated onto xdata. Dispatch is driven by
%   T.Properties.CustomProperties.Kind; on any evaluation error a matrix
%   of Inf is returned so calling optimizers can reject the point
%   instead of erroring out.
%
%   Inputs:
%     par   <-- row vector (or N x Np matrix for Kind==1, Simulink) of
%               parameter values to evaluate
%     T     <-- summary table from gsua_dataprep (or its
%               CustomProperties struct directly)
%     xdata <-- domain at which to evaluate/interpolate the model output
%
%   Outputs:
%     y <-- model output evaluated at xdata (Inf-filled if evaluation
%           failed)
%
%   See also GSUA_PARDEVAL, GSUA_EVAL, GSUA_DATAPREP.
try
    Table=Table2.Properties.CustomProperties;
catch
    Table=load('ATable.mat');
    Table=Table.Table2;
end
kind=Table.Kind;   
  
if kind==1
    [y,~]=sens_montecarlo(Table2,par,xdata,false,[0,size(par,1)]);
    if size(y,1)==1
        y=squeeze(y)';
    end
    return
end
ni=Table.NumVars;
fun=Table.Solver; 
out=Table.output;
try
switch kind
    case {2,6}
        y=fun(par);
    case 3
        y=deval(fun(par(:,1:ni),par(:,ni+1:end)),xdata,out);
    case 4
        y=gsua_intrp(fun(par(:,1:ni),par(:,ni+1:end)),xdata,out);
    case 5
        domain=Table.Domain;
        opt=Table.copt;
        y=gsua_intrp(fun(par,domain,opt),xdata,out,2);
end
catch ME
    warning(ME.message)
    y=inf(length(out),max(size(xdata))); 
end
end