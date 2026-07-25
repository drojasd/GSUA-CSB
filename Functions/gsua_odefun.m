function f = gsua_odefun(expr, vars, varargin)
%ODEFUNCTION   Convert a system (vector or matrix) of symbolic algebraic 
%   expressions to a MATLAB function handle suitable for ODE45, ODE15s,
%   and other ODE solvers.
%
%   f = odeFunction(expr, vars, options)
%   f = odeFunction(expr, vars, p1, p2, ..., options)
%
%   expr is a vector or matrix of symbolic expressions in the variables 
%   specified by the vector vars, which can consist of symbolic variables, 
%   univariate symbolic functions, or univariate symbolic function calls. 
% 
%   If expr contains symbolic parameters other than the variables given 
%   in vars, you must specify these parameters as additional arguments 
%   p1,p2,... .
%   These parameters can be symbolic variables, univariate or multivariate
%   symbolic functions, or univariate or multivariate symbolic function 
%   calls.
%
%   You can also use the name-value pair arguments acceptable by
%   MATLABFUNCTION, except that you cannot use 'vars' to control the 
%   order of the input variables. ODEFUNCTION passes these name-value
%   pair arguments to an internal MATLABFUNCTION call.
%
%   The output f is a function handle to a function f(t,Y,p1,p2,...).
%   The scalar input parameter t represents the 'time' variable. The
%   (univariate) symbolic functions or symbolic function calls listed
%   in vars depend on this variable.
%   The column vector Y represents the state variables (functions of 
%   'time') given by the input argument vars.
%   The arguments p1,p2,... represent the symbolic input parameters 
%   (any symbolic variables that do not appear in vars).
% 
%   Given numerical values v1,v2,... for the symbolic parameters 
%   p1,p2,..., the reduced function handle
%     F = @(t,Y) f(t,Y,v1,v2,...)
%   represents the function F(t,Y) = expr. The reduced function handle
%   F is suitable as the first argument to various numeric ODE solvers 
%   such as ODE45, ODE15s and so on.
%
%   Example:
%     Consider the following system of differential algebraic equations 
%     (DAEs):
%
%     >> syms x1(t) x2(t) a b r(t);
%        expr = [diff(x1(t),t) == a*x1(t) + b*x2(t)^2,... 
%                x1(t)^2 + x2(t)^2 == r(t)^2];
%        vars = [x1(t),x2(t)]
%
%     x1(t) and x2(t) are the state variables of the system. The system
%     also involves the parameters a and b, and the parameter function 
%     r(t). 
%     Extract the mass matrix of this DAE system:
%     
%     >>  [m,f] = massMatrixForm(expr,vars) 
%
%     Convert the expressions m and f to function handles depending on 
%     the variables x1, x2 and the parameters a, b, and r(t):
%
%     >> m = odeFunction(m,vars,a,b,r(t)) 
%     >> f = odeFunction(f,vars,a,b,r(t)) 
%   
%     Before using the numerical MATLAB ODE solvers, you must assign
%     numerical values to the parameters of the system.
%
%     >> a = -0.6;
%        b = -0.1;
%        r = @(t) cos(t)/(1+t^2);
%
%     Find the reduced function handles M and F.
%
%     >> M = @(t,Y) m(t,Y,a,b,r(t));
%     >> F = @(t,Y) f(t,Y,a,b,r(t));
%
%     Find consistent initial conditions for this DAE system:
%
%     >> t0 = 0;
%        y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
%        yp0= [a*y0(1) + b*y0(2)^2; 1.234];
%
%     Set the mass matrix and call ODE15S with these initial conditions.
%
%     >> opt = odeset('mass',M,'InitialSlope',yp0);
%     >> ode15s(F,[t0,1],y0,opt)
%         
% See also: MASSMATRIXFORM, ODE45, ODE15S.

%   Copyright 2014-2016 The MathWorks, Inc.

  narginchk(2, Inf);
  if isa(expr,'symfun')
    expr = formula(expr);
  end
   if isa(vars,'symfun')
    vars = formula(vars);
  end
%     
  isvect = isvector(expr);
  s = size(expr);
%   [expr, vars] = checkDAEInput(expr, vars);

  % Look for the position of the first string argument.
  % This is where the matlabFunction options should start:
  flagpos = find(cellfun(@(x) ischar(x) || isstring(x),varargin),1);
  if isempty(flagpos)
     params = varargin;
     matlabFunOptions = {};
  else
     params = varargin(1:flagpos-1);
     matlabFunOptions = varargin(flagpos:end);
  end
  % Do not allow the name 'vars' in the matlabFunction options:
  if any(cellfun(@(x) (ischar(x) || isstring(x)) && ...
     ~isempty(regexpi(x,'^(v|va|var|vars)$')),matlabFunOptions(1:2:end)))
     error(message('symbolic:odeFunction:UnexpectedVars'));
  end
  
  % check that the parameters are syms
  params = cellfun(@sym, params, 'UniformOutput', false);
 
  A = feval(symengine, 'daetools::odeFunction', expr, vars, params{:});
  A = num2cell(A);

  expr = A{1}.';
  % if expr is a matrix (e.g., the mass matrix),
  % then we need to restore the matrix form:
  if ~isvect
    expr = reshape(expr, s);
  end 
  
  if numel(params) == 1
      A{4} = A(4);
  else
      A{4} = num2cell(A{4});
  end

  varsAndParams = [{A{2}, A{3}} A{4}];
  f = matlabFunction(expr, 'Vars', varsAndParams, matlabFunOptions{:});
end