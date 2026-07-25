function [depth,idx]=gsua_depth(x,parallel)
%GSUA_DEPTH Generalized band depth of a set of curves, ranked deepest-first
%
%   [depth,idx] = gsua_depth(x)
%   [depth,idx] = gsua_depth(x,parallel)
%
%   Computes the generalized band depth of every curve in x relative to
%   the full set x (Lopez-Pintado & Romo style band depth), a measure of
%   how "central" each curve is among the others. Used by GSUA_SA in the
%   multi-objective sensitivity analysis branch to automatically pick a
%   representative reference output (the deepest curve) when no
%   experimental/nominal output ('ynom') is supplied.
%
%   Inputs:
%     x        <-- n x d x inputs array: n curves, d sample points per
%                  curve, inputs output signals (depth is averaged
%                  across signals)
%     parallel <-- (optional, logical) use parfor for the depth
%                  computation. Default: false
%
%   Outputs:
%     depth <-- n x 1 vector of average band depth per curve
%     idx   <-- indices of the curves sorted from deepest (most central,
%                idx(1)) to shallowest
%
%   See also GSUA_SA.
if nargin<2
    parallel=false;
end
if parallel
    parforArg=Inf;
else
    parforArg=0;
end

[n,d,inputs]=size(x);

depth=zeros(n,inputs);

for h=1:inputs
    PosiXX=zeros(n,d);
    parfor (i=1:n,parforArg)
        for j=1:d
            PosiXX(i,j)=sum(x(:,j,h)<x(i,j,h))+1;
        end
    end
    %disp(PosiXX)
    %disp(x);
    contador=n*(n+1)/2;
    parfor (i=1:n,parforArg)
        for j=1:d
            depth(i,h)=depth(i,h)+(PosiXX(i,j)-1)*(n-PosiXX(i,j))+n;
        end
    end

    depth(:,h)=depth(:,h)/contador;
end
depth=mean(depth,2);
[~,idx]=sort(depth,'descend');
end
        