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
%   Margin convention: the regulator GSUA_PE builds is mean((ydata - ydata*(m+1)).^2), which equals
%   (m*rms(ydata)).^2 per row, where m is the RAW margin you passed GSUA_PE ('margin',m). So with a
%   good-shape fit (correlation ~ 1) this cost falls below 1 exactly when the residual RMS is within
%   a fraction m of the data RMS (e.g. m = 0.3 -> "within ~30%"). Because the regulator does not
%   depend on the parameters, m only RESCALES the cost -- it does not move the optimum, but it sets
%   what "cost < 1" and the GSUA_NOISEFLOOR threshold mean. See GSUA_PE's 'margin' help.
%
%   See also GSUA_PE, GSUA_RCOSTF, GSUA_COSTFMULTI, GSUA_CORR2OMITNAN.
nvalid=sum(~isnan(ydata),2);
cost=(sum((ydata-yfunction).^2,2,'omitnan')./nvalid)./regulator;
% A row with a non-finite/non-positive regulator (an all-zero or constant data row, which GSUA_PE
% flags and sets to Inf) carries no usable scale and has an undefined correlation; exclude it
% entirely -- contribute 0 and drop it from the average -- instead of letting NaN poison the scalar
% cost (which otherwise makes fmincon report "objective undefined at initial point"). With every row
% finite (the usual case) keep is all-true and this is identical to the previous behaviour.
keep = isfinite(regulator) & regulator>0;
for i=1:inputs
    if ~keep(i)
        cost(i)=0;
        continue
    end
    cost(i)=((2-gsua_corr2omitnan(ydata(i,:),yfunction(i,:)))*cost(i))^alpha;
end
cos=sum(cost)/max(sum(keep),1);
end