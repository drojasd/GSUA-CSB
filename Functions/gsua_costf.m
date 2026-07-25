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
%   Inputs:
%     inputs    <-- number of output rows (signals) being compared
%     regulator <-- 1 x inputs normalization vector, typically the MSE
%                   between ydata and a margin-perturbed ydata (see
%                   gsua_pe)
%     len       <-- number of samples per output row
%     ydata     <-- inputs x len matrix of reference/experimental data
%     yfunction <-- inputs x len matrix of model output for the
%                   candidate parameter set
%     alpha     <-- exponent applied to each per-output cost term
%
%   Outputs:
%     cos <-- scalar cost (mean of the per-output penalized costs)
%
%   See also GSUA_PE, GSUA_RCOSTF, GSUA_COSTFMULTI.
cost=(sum((ydata-yfunction).^2,2)/len)./regulator;
for i=1:inputs
    cost(i)=((2-corr2(ydata(i,:),yfunction(i,:)))*cost(i))^alpha;
end
cos=sum(cost)/inputs;
end