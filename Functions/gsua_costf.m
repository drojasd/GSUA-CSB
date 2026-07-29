function cos=gsua_costf(inputs,regulator,len,ydata,yfunction,alpha)
%GSUA_COSTF Scalar correlation-penalized MSE cost function
%
%   cos = gsua_costf(inputs,regulator,len,ydata,yfunction,alpha)
%
%   Internal objective function used by GSUA_PE (for the 'ga',
%   'particle', 'psearch', 'surrogate', 'annealing' and 'fmincon'
%   solvers) to score a candidate parameter set. For each output row,
%   the mean squared error between ydata and yfunction is normalized by
%   regulator and penalized by a correlation factor (2-corr2), then
%   raised to alpha; the per-output costs are averaged into a single
%   scalar so any MATLAB scalar optimizer can be used.
%
%   NaN-tolerant per output row: a NaN in ydata (e.g. a gap in one
%   signal's real-world sampling, for a multi-output joint fit where
%   different signals have independently missing periods) is excluded
%   from that row's MSE and correlation terms rather than propagating
%   through and poisoning the whole scalar cost. The per-row MSE is
%   normalized by that row's own count of non-NaN samples, not the
%   fixed len, so a gappy row isn't silently deflated relative to a
%   fully-populated one. Identical to the pre-NaN-tolerant behavior when
%   ydata has no NaN (the per-row non-NaN count then equals len exactly).
%
%   Inputs:
%     inputs    <-- number of output rows (signals) being compared
%     regulator <-- 1 x inputs normalization vector, typically the MSE
%                   between ydata and a margin-perturbed ydata (see
%                   gsua_pe)
%     len       <-- number of samples per output row (retained for
%                   signature compatibility; the per-row MSE divisor is
%                   now each row's own non-NaN count, see above)
%     ydata     <-- inputs x len matrix of reference/experimental data,
%                   may contain NaN for missing samples
%     yfunction <-- inputs x len matrix of model output for the
%                   candidate parameter set
%     alpha     <-- exponent applied to each per-output cost term
%
%   Outputs:
%     cos <-- scalar cost (mean of the per-output penalized costs)
%
%   See also GSUA_PE, GSUA_RCOSTF, GSUA_COSTFMULTI, GSUA_CORR2OMITNAN.
nvalid=sum(~isnan(ydata),2);
cost=(sum((ydata-yfunction).^2,2,'omitnan')./nvalid)./regulator;
for i=1:inputs
    cost(i)=((2-gsua_corr2omitnan(ydata(i,:),yfunction(i,:)))*cost(i))^alpha;
end
cos=sum(cost)/inputs;
end