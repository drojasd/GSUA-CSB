function sol = orientationFixtureModel(pars, domain, config)
%ORIENTATIONFIXTUREMODEL Test fixture: 3-state linear decay in selectable sol orientation.
%
%   sol = orientationFixtureModel(pars, domain, config)
%
%   Integrates the decoupled linear system  dy_j/dt = -k_j*y_j  (analytic
%   solution y_j(t) = y0_j*exp(-k_j*t)) and returns the ODE solution struct
%   in whatever orientation config.form requests, to exercise the
%   orientation normalization at the top of gsua_intrp:
%
%     'canonical'        -> sol.x = 1xN, sol.y = nStates x N, sol.solver = ''
%     'notsolver'        -> native nStates x N but sol.solver = 'ode45'
%     'ty'               -> the [t,y]=ode45 mistake: sol.x = Nx1, sol.y = N x nStates,
%                           no sol.solver field
%     'structTransposed' -> the Yang mistake: sol.y transposed to N x nStates
%                           with sol.solver = 'ode45'
%
%   pars = [y0(1:3); k(1:3)].  nStates = 3, deliberately far below the number
%   of ode45 time steps so orientation is never ambiguous.
y0 = pars(1:3);
k  = pars(4:6);

form = 'canonical';
if nargin >= 3 && isstruct(config) && isfield(config,'form')
    form = config.form;
end

[t,y] = ode45(@(tt,yy) -k(:).*yy, domain, y0(:));   % t: Nx1 , y: N x 3

switch form
    case 'canonical'
        sol.x = t.'; sol.y = y.'; sol.solver = '';
    case 'notsolver'
        sol.x = t.'; sol.y = y.'; sol.solver = 'ode45';
    case 'ty'
        sol.x = t;   sol.y = y;                       % Nx1 and N x nStates
    case 'structTransposed'
        sol.x = t.'; sol.y = y;   sol.solver = 'ode45';   % y is N x nStates
    otherwise
        error('orientationFixtureModel:BadForm','Unknown form "%s".', form);
end
end
