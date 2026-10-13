function r = gsua_corr2(a,b)
%GSUA_CORR2 2-D (Pearson) correlation coefficient without Image Processing Toolbox.
%
%   r = gsua_corr2(a,b)
%
%   Returns the correlation coefficient between two equally sized arrays a
%   and b, computed over all their elements. This reproduces the Image
%   Processing Toolbox function CORR2 exactly -- including returning NaN
%   (0/0) when either input is constant -- so the toolbox's cost and
%   sensitivity functions no longer depend on that product.
%
%   Inputs:
%     a,b <-- equal-size vectors or matrices to correlate
%
%   Outputs:
%     r <-- correlation coefficient in [-1,1], or NaN if an input is constant
%
%   See also GSUA_CORR2OMITNAN.
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
r = sum(a.*b) / sqrt(sum(a.^2) * sum(b.^2));
end
