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
if extra==1
    if ~isa(sol,'struct')  
      % Try  DEVAL(XINT,SOL)
      temp = struct();
      sol  = xint;
      xint = temp;
    else
        try
          t = sol.x;    
          y = sol.y; 
          if isfield(sol,'solver')
              solver=sol.solver;
          else
              solver='';
          end
          reps=length(idx);
          if any(strcmp(solver,{'dde23',''}))
              a=zeros(1,size(y,1),size(y,2));
              a(1,:,:)=y;
              y=a;
          end

        catch
          error(message('MATLAB:deval:SolNotFromDiffEqSolver', inputname( 1 )));
        end

    end


    sz=size(y);
    if nargin < 3
      idx = 1:sz(2);  % return all solution components
    else 
      if any(idx < 0) || any(idx > size(y,2))
        error(message('MATLAB:deval:IDXInvalidSolComp', inputname( 3 )));
      end  
    end  

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