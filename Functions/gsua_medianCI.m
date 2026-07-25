function [CI_lower, CI_upper] = gsua_medianCI(X, alpha)
% BINOMIALMEDIANCI  Compute a nonparametric (distribution-free) confidence
% interval for the median of each row of X, using the order-statistics and 
% binomial approach.
%
%   X is an M-by-N matrix of data. Each row is a separate variable
%   (N samples per variable).
%
%   ALPHA is the significance level (e.g., 0.05 for a 95% confidence interval).
%   If ALPHA is omitted, it defaults to 0.05.
%
%   [CI_LOWER, CI_UPPER] = BINOMIALMEDIANCI(...) returns two M-by-1 vectors 
%   containing the lower and upper confidence limits for each row's median.
%
%   This method does NOT assume normality; it uses the distribution-free 
%   construction:
%       CI = [ X_(r) ,  X_(s) ]
%   where X_(k) denotes the k-th order statistic, and r and s are chosen
%   by inverting a Binomial(n, 1/2) distribution to achieve coverage
%   approximately 1 - ALPHA.
%
%   See also GSUA_IA, GSUA_DIA.

    if nargin < 2
        alpha = 0.05;  % default for a 95% CI
    end
    
    [M, N] = size(X);
    
    % Preallocate outputs
    CI_lower = NaN(M, 1);
    CI_upper = NaN(M, 1);
    
    % We want indices r and s so that:
    %   P(X_(r) <= true median <= X_(s)) ≈ 1 - ALPHA
    % The usual approach:
    %   r = smallest integer such that sum_{k=0}^{r-1} [C(n,k)*(1/2)^n] >= ALPHA/2
    %   s = largest integer  such that sum_{k=0}^{s-1} [C(n,k)*(1/2)^n] <= 1 - ALPHA/2
    
    % We'll figure out r and s *once* based on N (assuming all rows have the same sample size).
    
    % We can use binocdf(k, N, 0.5) = sum_{i=0}^k C(N,i)*(1/2)^N.
    
    % --- Find r ---
    r = 1;
    while r <= N && binocdf(r-1, N, 0.5) < alpha/2
        r = r + 1;
    end
    % Safety clamp in case the loop never stops before exceeding N
    if r > N
        r = N;
    end
    
    % --- Find s ---
    s = N;
    while s >= 1 && binocdf(s-1, N, 0.5) > 1 - alpha/2
        s = s - 1;
    end
    % Safety clamp
    if s < 1
        s = 1;
    end
    
    % For each row (variable), sort and pick out X_(r) and X_(s)
    for i = 1:M
        rowData = sort(X(i, :));
        CI_lower(i) = rowData(r);
        CI_upper(i) = rowData(s);
    end
end
