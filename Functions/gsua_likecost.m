function cost=gsua_likecost(ydata,ynom,margin)
%GSUA_LIKECOST Gaussian negative-log-likelihood cost relative to a margin baseline
%
%   cost = gsua_likecost(ydata,ynom,margin)
%
%   Standalone likelihood-style objective function: computes the
%   Gaussian negative log-likelihood of each run in ydata against ynom
%   (assuming a standard deviation proportional to ynom*margin), then
%   subtracts the negative log-likelihood of the margin-perturbed
%   baseline itself. The result is a relative cost that is zero when a
%   run matches ynom as closely as the margin baseline does, and grows
%   as the run deviates further. This mirrors the margin<0 objective
%   built inline in GSUA_PE and can be used as a standalone objective
%   function for custom optimization setups.
%
%   Inputs:
%     ydata  <-- reps x len x 1 array, one row per simulation run
%     ynom   <-- 1 x len reference/nominal output
%     margin <-- relative standard deviation used for both the
%                likelihood weighting and the baseline (absolute value
%                is used internally)
%
%   Outputs:
%     cost <-- 1 x reps row vector of relative negative log-likelihood
%              costs, one per run
%
%   See also GSUA_PE, GSUA_LIKELIHOOD.

margin=1+abs(margin);
ymargin=ynom*margin;
desv = abs((ymargin).^2);
[reps,~,~]=size(ydata);
cost = zeros(1,reps);
costmargin = sum(log(2*pi*desv) +(ymargin-ynom).^2./desv,'all','omitnan')/2;

for j=1:reps
    cost(j)=sum(log(2*pi*desv) +(ydata(j,:,1)-ynom).^2./desv,'all','omitnan')/2;
end
cost=cost-costmargin;
end