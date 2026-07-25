function cos=gsua_costfMulti(ydata,yfunction,margin,parallel,alpha)
%GSUA_COSTFMULTI Vectorized correlation-penalized MSE cost over many simulations
%
%   cos = gsua_costfMulti(ydata,yfunction)
%   cos = gsua_costfMulti(ydata,yfunction,margin,parallel,alpha)
%
%   Batched version of GSUA_COSTF: scores every simulation ("run") in a
%   3-D array of model outputs against a single reference yfunction in
%   one call. Used by GSUA_SA (multi-objective sensitivity analysis) and
%   GSUA_CSB to reduce a whole population of Monte-Carlo/design-matrix
%   simulations to one scalar cost per run.
%
%   Inputs:
%     ydata     <-- reps x len x inputs array, one page of len-sample
%                   output per input signal, one row per simulation run
%                   (typically the output of gsua_pardeval)
%     yfunction <-- inputs x len reference/nominal output (or a
%                   1 x len x inputs array reshaped internally)
%     margin    <-- (optional) relative perturbation used to build the
%                   internal per-signal regulator (values <=1 are
%                   treated as 1+abs(margin)). Default: 0.1
%     parallel  <-- (optional, logical) evaluate runs with parfor.
%                   Default: false
%     alpha     <-- (optional) exponent applied to each per-signal cost
%                   term. Default: 1
%
%   Outputs:
%     cos <-- 1 x reps row vector of scalar costs, one per run (Inf for
%             runs where the cost computation failed)
%
%   See also GSUA_COSTF, GSUA_SA, GSUA_CSB, GSUA_PARDEVAL.
if nargin<3
    margin=0.1;
    parallel=false;
end
if nargin<5
    alpha=1;
end
if parallel
    parforArg=Inf;
else
    parforArg=0;
end
if margin<=1
    margin=1+abs(margin);
end

[reps,len,inputs]=size(ydata);
if size(yfunction,3)>1
    yfunction=squeeze(yfunction)';
end
regulator=sum((yfunction-yfunction*margin).^2,2)/len;
cos=zeros(1,reps);

parfor (j=1:reps,parforArg) 
%for j=1:reps
    cost=zeros(1,inputs);
    try
        for i=1:inputs
            icorre=(2-corr2(ydata(j,:,i),yfunction(i,:)));
            if isnan(icorre)
                icorre=1;
            end
            cost(i)=(icorre*(sum((ydata(j,:,i)-yfunction(i,:)).^2,2)/len)/regulator(i))^alpha;
        end
        cos(j)=sum(cost)/inputs;
    catch
        cos(j)=inf;
    end
    
%     if isnan(cos(j))
%         disp('fail')
%     end
end
end