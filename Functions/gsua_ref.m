function optim=gsua_ref(est,table,xdata,ydata,exe)
%GSUA_REF Iteratively re-estimate a parameter set while nudging ranges away from their bounds
%
%   T = gsua_ref(est,T,xdata,ydata,exe)
%
%   Repeatedly re-runs GSUA_PE with fmincon starting from est: whenever
%   a (non-excluded) parameter's current estimate sits within 10% of
%   its range bound, that bound is stretched (widened by 50% on the
%   binding side) and the parameter set is re-estimated. This continues
%   until every non-excluded parameter estimate settles comfortably
%   inside its range. Useful as a quick range-adequacy check/refinement
%   after an initial GSUA_PE call, when you suspect the parameter
%   bounds are actively constraining the fit rather than the model
%   itself.
%
%   Inputs:
%     est   <-- initial parameter estimate (row vector, one value per
%               table row)
%     table <-- summary table from gsua_dataprep (or a prior gsua_pe
%               result); table.Range is expanded in place as needed
%     xdata <-- domain at which the model is evaluated
%     ydata <-- experimental/reference output at xdata
%     exe   <-- indices of parameters to exclude from the bound check
%               (kept fixed at their current range)
%
%   Outputs:
%     optim <-- updated summary table after convergence, with Range
%               expanded as needed and a final GSUA_PE estimation column
%
%   See also GSUA_PE, GSUA_OATR, GSUA_OATR2.
ex=true(1,size(table,1));
ex(exe)=false;
solver='fmincon';
opt=optimoptions('fmincon','UseParallel',false,'Display','iter','MaxFunctionEvaluations',5000);
flag=true;
while flag
    range=table.Range;
    if any(any([range(ex,2)<est(ex)*1.1,range(ex,1)>est(ex)*0.9]))
        A=find(range(:,2)<est*1.1);
        B=find(range(:,1)>est*0.9);
        A=setdiff(A,exe);
        B=setdiff(B,exe);
        table.Range(A,2)=est(A)*1.5;
        table.Range(B,1)=est(B)*0.5;
    else
        break
    end

    [table,~]=gsua_pe(table,xdata,ydata,'solver',solver,'N',1,'opt',opt,'ipoint',est');
    est=table.Estfmincon(:,1);
end
optim=table;
end
