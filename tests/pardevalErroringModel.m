function sol = pardevalErroringModel(pars, domain, ~)
%PARDEVALERRORINGMODEL Test fixture: a Kind-5 model that always errors.
%
%   Used by gsuaPardevalFailWarningTest to check that gsua_pardeval reports failed
%   simulations instead of silently leaving them as zeros.
error('pardevalErroringModel:boom', 'deliberate failure for testing'); %#ok<NASGU>
end
