function [C,E]=GeometricMedian(X,Co,opt,W)
% This function takes as input a d-dimensional (d>1) dataset composed on N
% samples and computes its geometric median. Geometric median of a 
% discrete set of sample points in Euclidean space is defined as the 
% point minimizing the sum of ABSOLUTE distances to the sample 
% points [1]. By contrast, the arithmetic mean is the point minimizing the
% sum of SQUARED distances. In comparison to the arithmetic mean, geometric
% median provides a more robust estimate of the central tendency of the 
% data (i.e., mean can be affected by an arbitrary amount even if there is 
% a single outlier, and therefore has a break-down point of 0%, where as 
% median has a break-down point of 50%). Since closed-form formula for 
% computing the (geometric) median does not exist, this statistic must 
% be estimated using iterative optimization methods. Here, the geometric 
% median is computed using Weiszfeld's algorithm [1].
%
% INPUT:
%   - X     : N-by-d array of d-dimensional data points/vectors, where N is
%             the total number of samples.
%   - Co    : optional input argument specifying the starting point for the
%             optimization. Co=(W'*X)/sum(W) is the default setting; see 
%             definition of W below.
%   - opt   : optional input argument specifying converge criteria;
%             opt=[Nmax tol], where Nmax is maximum number of iterations
%             and tol is maximum change in position of the median 
%             between two successive iterations. opt=[50 1E-6] is the 
%             default setting. Optimization terminates when either one
%             of the above criteria is met.
%   - W     : optional input argument. W is a N-by-1 vector of (positive) 
%             weights assigned to the points in X. W=ones(N,1)/N is the 
%             default setting.
%
% OUTPUT:
%   - C     : 1-by-d vector specifying geometric median of X.     
%   - E     : 1-by-(K+1) vector containing values of the total (weighted) 
%             absolute distance from X to C_k where C_k is the estimate of
%             C at iteration k; K is the total number of iterations.
%             E(1) corresponds to initialization.
%
% REFERENCES:
% [1] http://en.wikipedia.org/wiki/Geometric_median
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<3 || isempty(opt), opt=[50 1E-6]; end
opt=abs(opt);

if ~isnumeric(X) || ~ismatrix(X)
    error('1st input argument (X) must be a 2D array, with observations along the rows')
end

if numel(opt)~=2 || ~isnumeric(opt)
    error('Converge criteria must be specified as a 1-by-2 array; [Nmax tol]. See function description for more info.');
end

E=[];
if isempty(X), C=[]; return; end
d=size(X,2);

if nargin<4 || isempty(W)
    W=ones(size(X,1),1); 
elseif numel(W)~=size(X,1) || sum(W<0)>1
    error('Invalid format for 4th input argument (W)')
end
W=abs(W(:));
W=W/sum(W);

if nargin<2 || isempty(Co)
    Co=W'*X;
end
if numel(Co)~=d
    error('Dimensionality of the starting point does not match dimensionality of the data')
end
Co=Co(:)';
if d==1, Co=median(X); end


% Sum of distances
if nargout>1
    E=W'*sqrt(sum(bsxfun(@minus,X,Co).^2,2));
end

% Compute geometric median
C=Co; dC=Inf; opt(2)=max(opt(2).^2,1E-16); 
a=1E-1;
n=1;
while n<=opt(1) && dC>opt(2)
    
    n=n+1;
    
    w=sqrt(sum(bsxfun(@minus,X,C).^2,2))./W;
    
    if nargout>1, E(1,n)=sum(w); end %#ok<*AGROW>
    
    w=1./(w+a); % a is added for 2 reasons: 1) to avoid potential division by 0, and 2) to help overcome local minima when C is close to one of the sample points 
    
    Cn=sum(bsxfun(@times,X,w),1)/sum(w);
    
    dC=sum((C-Cn).^2);
    C=Cn;
    a=max(a/10,eps); % relax regularization parameter a

    %fprintf('%3u   %.3E \n',n,sqrt(dC/opt(2)))
end

