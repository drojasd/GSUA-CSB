function [Table,solver] = gsua_dpmat(odes,vars,domain,modelName,varargin)
% Consolidation of data and symbolic matlab environment for (c) GSUA-RC Toolbox
%
% [T,solver]=gsua_dpmat(odes,vars,domain,model_name)
% Parameters:
% odes       <-- array with the system of differential equations
% vars       <-- array of state variables for the system of equations
% domain     <-- array with the domain of the system of equations
% model_name <-- name for the output file of the system
% Outputs:
% T      <-- table array with all required system information
% solver <-- handle function for ode solvers
%
% To explore additional configuration, open the respective user guide:
% open gsua_matuserguide <--- symbolic matlab userguide

%input parser scheme
p=inputParser;
checkOdes=@(x) isa(x,'symfun')||isa(x,'sym');
checkDomain=@(x) isnumeric(x) && size(x,2)<=3;

defaultrMethod='range';
validrMethod={'percent','range','normal', 'std'};
checkrMethod = @(x) any(validatestring(x,validrMethod));

defaultRange=[];
checkRange = @(x) (isnumeric(x) || isempty(x));
checkModelName=@(x) (isempty(x) || ischar(x));
defaultNominals=[];
defaultOutput=[];
dfsolver='ode45';
validSolver={'ode45','ode4'};
checkSolver=@(x) any(validatestring(x,validSolver));
dfOpt=[];

addRequired(p,'odes',checkOdes);
addRequired(p,'vars',checkOdes);
addRequired(p,'domain',checkDomain);
addRequired(p,'modelName',checkModelName);
addParameter(p,'range',defaultRange,checkRange);
addParameter(p,'rMethod',defaultrMethod,checkrMethod);
addParameter(p,'nominal',defaultNominals,@isnumeric);
addParameter(p,'output',defaultOutput,checkRange);
addParameter(p,'opt',dfOpt);
addParameter(p,'solver',dfsolver,checkSolver);


parse(p,odes,vars,domain,modelName,varargin{:})
odes=p.Results.odes;
vars=p.Results.vars;
domain=p.Results.domain;
Range=p.Results.range;
rMethod=p.Results.rMethod;
nominal=p.Results.nominal;
model_name=p.Results.modelName;
output=p.Results.output;
opt=p.Results.opt;
sname=p.Results.solver;
pODEs = symvar(odes);
pODEvars = symvar(vars);
extraPars = setdiff(pODEs, pODEvars);
tPars = [vars,extraPars];
tPars=char(tPars);
tPars=split(tPars(10:end-3),', ');
Np=max(size(tPars));

% Output table construction
if isempty(Range)
    Range=zeros(Np,2);
    special=false;
else
    special=true;
    if max(size(Range)) ~= Np || min(size(Range)) ~= 2
        disp('Please provide ranges as a matrix of Nfactors*2')
        return
    else
        if size(Range,2)==Np
            Range=Range';
        end
    end
end
fixvars=[];
%fixed=[];
jump=true;

switch rMethod
    case {'range','normal'}
       % Table.Range=Range;
    case 'percent'
        Range=[Range(:,1).*(1-Range(:,2)/100), Range(:,1).*(1+Range(:,2)/100)];
    case 'std'
        Range=[Range(:,1)-Range(:,2), Range(:,1)+Range(:,2)];
    otherwise
        disp('Available rMethods are: ')
        dis('range, percent, normal, std')
        return
end
Table = table(Range);
if size(domain,2)==1
    NumVars=0;
    limit=0;
    jump=false;
    Table.Properties.RowNames=tPars;
else
    NumVars=max(size(split(char(vars),' ')));
    RNames=split(tPars(1:NumVars),'(');
    Table.Properties.RowNames = [strcat(RNames(:,1),'0');tPars(NumVars+1:end)];
end

if ~isempty(nominal)
    if max(size(nominal)) == Np && min(size(nominal)) == 1
        if max(size(nominal))==size(nominal,2)
            Table.Nominal=nominal';
        else
            Table.Nominal=nominal;
        end
    else
        disp('Wrong number of nominal values')
    end  
else
    switch rMethod
        case 'normal'
            Table.Nominal=Table.Range(:,1);
        otherwise
            Table.Nominal=(Table.Range(:,1)+Table.Range(:,2))./2;
    end
end

if special && jump
    fixvars=Table.Range(1:NumVars,1)'==Table.Range(1:NumVars,2)';
    limit=sum(fixvars);
    NumVars=NumVars-limit;
    
    if limit~=sum(fixvars(1:limit))
        disp('All fixed variables must be at begining!')
        return
    end
end
disp('Introduce ranges in the following order: ')
[vars,extraPars]
if isempty(extraPars)
    extraPars=vars;
end
if special
    fixed=Table.Range(NumVars+limit+1:end,1)'==Table.Range(NumVars+limit+1:end,2)';
    helper=false(1,NumVars);
    if any(fixed)
        odes=subs(odes,extraPars(fixed),Table.Range([helper,fixed],1)');
        pODEs = symvar(odes);%new parameters
        extraPars = setdiff(pODEs, pODEvars);%new parameters
    end
else
    fixed=false(1,length(vars));
end

%Function handle crafting
if size(domain,2)==1
    kind=2;
    fun=matlabFunction(odes,'file',genvarname(model_name),'vars',vars(~fixed));
    solver=@(pars) call_on(fun,pars);
    if isempty(output)
        output=1;
    end
    NumVars=size(output,2);
    RNames=cellstr(char(string(1:NumVars)'));
else
    [Mass,F] = massMatrixForm(odes,vars);
    fun=Mass\F;
    if strcmp(sname,'ode45')
        odefun=odeFunction(fun,vars,extraPars,'File',genvarname(model_name));
        solver= @(init,pars) ode45(@(t,Y) odefun(t,Y,pars),domain,[Table.Range(fixvars,1)',init],opt);
        kind=3;
    else
        odefun=gsua_odefun(fun,vars,extraPars,'File',genvarname(model_name));
        solver= @(init,pars) ode4(@(t,Y) odefun(t,Y,pars),domain,Table.Range(fixvars,1)',init,opt);
        kind=4;
    end
end

try
    Table=addprop(Table,{'Solver' 'NumVars' 'output' 'Domain' 'rMethod' 'Kind' 'Fixed' 'Vars','Sname'},{'table' 'table' 'table' 'table','table' 'table' 'table', 'table','table'});
    Table.Properties.CustomProperties.NumVars=NumVars;
    Table.Properties.CustomProperties.Domain=domain;
    Table.Properties.CustomProperties.rMethod=rMethod;
    Table.Properties.CustomProperties.Solver=solver;
    Table.Properties.CustomProperties.Kind=kind;
    Table.Properties.CustomProperties.Vars=RNames(:,1)';
    Table.Properties.CustomProperties.Fixed=[];
    
    if isempty(output)
        Table.Properties.CustomProperties.output=1:Table.Properties.CustomProperties.NumVars;
    else
        Table.Properties.CustomProperties.output=output;
    end

catch
    disp('Your Matlab version does not support custom properties on tables!')
    disp('Creating auxiliar table (ATable) as a file, please do not erase ATable')
    Table2=table(NumVars);
    Table2.Domain=domain;
    Table2.rMethod=rMethod;
    Table2.Solver=solver;
    Table2.Kind=kind;
    Table2.Vars=RNames(:,1)';
    Table2.Fixed=false(1,size(Table,1)-sum([fixvars,fixed]));
     
     if isempty(output)
        Table2.output=1:Table2.NumVars;
    else
        Table2.output=output;
     end
     save('ATable','Table2');
end
Table(logical([fixvars,fixed]),:)=[];
end

function r = call_on(func, array)
  t = num2cell(array,1);
  r = func(t{:});
end

%Special extra functions
function r = ode4(odefun,tspan,fix,y0,opt)
%ODE4  Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.   
%
%   Example 
%         tspan = 0:0.1:20;[]
%         y = ode4(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%
reps=size(y0,1);
y0=[repmat(fix,reps,1),y0];
tspan=tspan(1):tspan(3):tspan(2);
h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  
% y0=y0';
try
 feval(odefun,tspan(1),y0);
catch E
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',E];
  error(msg);  
end 
if ~isempty(opt)
    neg=opt.NonNegative;
else
    neg=[];
end
if isempty(neg)
    special=false;
else
    special=true;
end
neq = size(y0);
N = length(tspan);
Y = zeros(neq(1),neq(2),N);
F = zeros(neq(1),neq(2),4);

Y(:,:,1) = y0;
for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1); %% Vectorization have a failure
  yi = Y(:,:,i-1);  
  F(:,:,1) = reshape(feval(odefun,ti,yi),neq(1),neq(2));
  F(:,:,2) = reshape(feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,:,1)),neq(1),neq(2));
  F(:,:,3) = reshape(feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,:,2)),neq(1),neq(2));
  F(:,:,4) = reshape(feval(odefun,ti+hi,yi+hi*F(:,:,3)),neq(1),neq(2));
  Y(:,:,i) = yi + (hi/6)*(F(:,:,1) + 2*F(:,:,2) + 2*F(:,:,3) + (F(:,:,4)));
  if special && any(any(Y(:,neg,i)<0))
      y_temp=Y(:,neg,i);
      y_temp(y_temp<0)=0;
      Y(:,neg,i)=y_temp;
  end
end
%Y = Y.';
r=struct;
r.solver='ode4';
r.x=tspan;
r.y=Y;

end