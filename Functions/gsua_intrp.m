function Sxint = gsua_intrp(sol,xint,idx,extra)
%GSUA_INTRP Interpolate ODE solver output at requested time points
%
%   Sxint = gsua_intrp(sol,xint,idx)
%   Sxint = gsua_intrp(sol,xint,idx,extra)
%
%   Internal helper (DEVAL-style) used to sample the time-domain output
%   of an ODE solution structure (from ode45 or the bundled ODE4
%   integrator) at the requested time points xint, selecting the state
%   components given by idx. Used by GSUA_DEVAL and GSUA_PARDEVAL for
%   Kind==4/5 (custom fixed-step or user-defined-with-domain models).
%
%   Inputs:
%     sol   <-- ODE solution struct with fields x (time) and y (states,
%               possibly with a leading simulation-run dimension)
%     xint  <-- time points to interpolate at
%     idx   <-- indices of the state/output components to extract
%     extra <-- (optional) 1 to use the struct/sol-array form (default),
%               2 to use the simpler per-index interp1 form
%
%   Outputs:
%     Sxint <-- interpolated output at xint for the requested components
%
%   See also GSUA_DEVAL, GSUA_PARDEVAL.
if nargin <4
    extra=1;
end

% DEVAL-style argument order gsua_intrp(xint,sol): accept it by swapping, so the
% solution struct is always in sol before normalization below.
if ~isa(sol,'struct') && isa(xint,'struct')
    [sol, xint] = deal(xint, sol);
end

% Normalize a user-supplied ODE solution to the toolbox's canonical
% orientation, so that user-defined models are robust to how ode45 was
% called. ode45 has two output forms with opposite orientations:
%   sol = ode45(...)      -> sol.x is 1xN, sol.y is nStates x N
%   [t,y] = ode45(...)    -> t is Nx1,  y is N x nStates
% Either one (transposed or not, with or without sol.solver set) is accepted
% here and coerced to sol.x = 1xN, sol.y = nStates x N, sol.solver = ''.
% This keys off numel(sol.x): whichever dimension of a 2-D sol.y matches the
% number of time points is the time axis. Symbolic models (Kind 3) never
% reach gsua_intrp (they use MATLAB's deval), so only user-defined function
% models (Kind 4/5) are affected. A square output (N == nStates) is left
% untouched because its orientation is genuinely ambiguous.
if isa(sol,'struct') && isfield(sol,'x') && isfield(sol,'y') && ismatrix(sol.y)
    n = numel(sol.x);
    if n > 0
        sol.x = reshape(sol.x,1,n);
        yy = sol.y;
        if isvector(yy)
            sol.y = reshape(yy,1,n);
        elseif size(yy,1) == n && size(yy,2) == n
            % ambiguous square output: leave as-is
        elseif size(yy,2) == n
            % already nStates x N
        elseif size(yy,1) == n
            sol.y = yy.';               % N x nStates -> nStates x N
        else
            error('GSUA:gsua_intrp:BadShape', ...
                ['Model output sol.y is %s but sol.x has %d time points; ' ...
                 'sol.y must be nStates-by-N or N-by-nStates so that one ' ...
                 'dimension matches numel(sol.x).'], mat2str(size(yy)), n);
        end
    end
    sol.solver = '';
end

if extra==1
    try
          t = sol.x;    
          y = sol.y; 
          if isfield(sol,'solver')
              solver=sol.solver;
          else
              solver='';
          end
          if any(strcmp(solver,{'dde23',''}))
              a=zeros(1,size(y,1),size(y,2));
              a(1,:,:)=y;
              y=a;
          end

    catch
        error(message('MATLAB:deval:SolNotFromDiffEqSolver', inputname( 1 )));
    end


    sz=size(y);
    if nargin < 3
      idx = 1:sz(2);  % return all solution components
    else 
      if any(idx < 0) || any(idx > size(y,2))
        error(message('MATLAB:deval:IDXInvalidSolComp', inputname( 3 )));
      end  
    end  

    reps=length(idx);   % after idx has its default, so the 2-argument call works
    if reps==1
        if sz(1)>1
            Sxint=interp1(t,squeeze(y(:,idx,:))',xint)';
        else
            Sxint=interp1(t,squeeze(y(:,idx,:))',xint);
        end
    else
        Sxint=zeros(sz(1),length(xint),reps);
        for i=1:reps
            Sxint(:,:,i)=interp1(t,squeeze(y(:,idx(i),:))',xint)';
        end
    end
else
    t = sol.x;
    y = sol.y;
    reps=length(idx);
    Sxint=zeros(reps,length(xint));
    for i=1:reps
        Sxint(i,:)=interp1(t,y(idx(i),:),xint);
    end
end