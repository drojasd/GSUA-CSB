function cos=gsua_rcostf(ydata,yfunction,margin,alpha)
%GSUA_RCOSTF Scalar correlation-penalized MSE cost with a self-derived regulator
%
%   cos = gsua_rcostf(ydata,yfunction)
%   cos = gsua_rcostf(ydata,yfunction,margin,alpha)
%
%   Standalone alternative to GSUA_COSTF for use as a scalar objective
%   function with MATLAB optimizers (fmincon, ga, particleswarm, ...).
%   Unlike GSUA_COSTF, the normalization ("regulator") is derived
%   internally from ydata itself (the MSE between ydata and a
%   margin-perturbed copy of ydata), so no external regulator needs to
%   be precomputed by the caller.
%
%   Inputs:
%     ydata     <-- inputs x len matrix of reference/experimental data
%     yfunction <-- inputs x len matrix of model output for the
%                   candidate parameter set
%     margin    <-- (optional) relative perturbation used to build the
%                   internal regulator (values <1 are treated as
%                   1+abs(margin)). Default: 1.1
%     alpha     <-- (optional) exponent applied to each per-output cost
%                   term. Default: 2
%
%   Outputs:
%     cos <-- scalar cost (mean of the per-output penalized costs)
%
%   See also GSUA_COSTF, GSUA_COSTFMULTI, GSUA_PE.
if nargin<3
    margin=1.1;
end
if nargin<4
    alpha=2;
end

margin=abs(margin);

if margin<1
    margin=margin+1;
end

[inputs,lon] = size(ydata);
regulator=sum((ydata-ydata*margin).^2,2)/lon;
cost=(sum((ydata-yfunction).^2,2)/lon)./regulator;


for i=1:inputs
    factor = 2-corr2(ydata(i,:),yfunction(i,:));
    if isnan(factor)
        factor = 1;
    end
    cost(i)=(factor*cost(i))^alpha;
end
cos=sum(cost)/inputs;
end