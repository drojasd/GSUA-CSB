function [Table,solver] = gsua_userdefined(func,Range,varargin)
%Summary of the function is still under development

%input parser scheme
p=inputParser;


checkDomain=@(x) isnumeric(x) && size(x,2)<=3;
defaultDomain=1;
defaultrMethod='range';
validrMethod={'percent','range','normal', 'std'};
checkrMethod = @(x) any(validatestring(x,validrMethod));
checkRange = @(x) (isnumeric(x) || isempty(x));
defaultNominals=[];
defaultOutput=1;
defaultNames={};
defaultONames={'out'};
dfOpt=[];
defaultVector=false;

addRequired(p,'func');
addRequired(p,'Range',checkRange);
addParameter(p,'names',defaultNames,@iscellstr);
addParameter(p,'domain',defaultDomain,checkDomain);
addParameter(p,'rMethod',defaultrMethod,checkrMethod);
addParameter(p,'nominal',defaultNominals,@isnumeric);
addParameter(p,'output',defaultOutput,checkRange);
addParameter(p,'out_names',defaultONames,@iscellstr);
addParameter(p,'vectorized',defaultVector,@islogical);
addParameter(p,'opt',dfOpt);

parse(p,func,Range,varargin{:})
func=p.Results.func;
Range=p.Results.Range;


names=p.Results.names;
domain=p.Results.domain;
rMethod=p.Results.rMethod;
nominal=p.Results.nominal;
output=p.Results.output;
out_names=p.Results.out_names;
vectorized=p.Results.vectorized;
opt=p.Results.opt;

%Table crafting

switch rMethod%Craft table.ranges
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
Np=size(Table,1);

if ~isempty(nominal)%craft table.nominal
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
if ~isempty(names)
    Table.Properties.RowNames = names;
else
    Table.Properties.RowNames = split(num2str(1:Np));
end

fixed=Table.Range(:,1)==Table.Range(:,2);%Identify fixed parameters


%Function handle crafting
user_func=str2func(func);
inputs=abs(nargin(user_func));
kind=5;
if sum(fixed)>0
input_control=zeros(1,Np);
input_control(fixed)=Table.Range(fixed,1);
    switch inputs
        case 1
            solver=@(pars) user_func(fixing(pars,input_control,fixed)');
            if vectorized
                kind=2;
            else
                kind=6;
            end
        case 2
            solver=@(pars) user_func(fixing(pars,input_control,fixed)',domain);
        otherwise
            solver=@(pars) user_func(fixing(pars,input_control,fixed)',domain,opt);
    end
else
    switch inputs
    case 1
        solver=@(pars) user_func(pars');
        if vectorized
                kind=2;
            else
                kind=6;
            end
    case 2
        solver=@(pars) user_func(pars',domain);
    otherwise
        solver=@(pars) user_func(pars',domain,opt);
    end
    
end

%Table custom properties
try
    Table=addprop(Table,{'Solver' 'NumVars' 'output' 'Domain' 'rMethod' 'Kind' 'Fixed' 'Vars',},{'table' 'table' 'table' 'table','table' 'table' 'table', 'table'});
    Table.Properties.CustomProperties.NumVars=length(out_names);
    Table.Properties.CustomProperties.Domain=domain;
    Table.Properties.CustomProperties.rMethod=rMethod;
    Table.Properties.CustomProperties.Solver=solver;
    Table.Properties.CustomProperties.Kind=kind;
    Table.Properties.CustomProperties.Vars=out_names;
    Table.Properties.CustomProperties.Fixed=[];
    
    if isempty(output)
        Table.Properties.CustomProperties.output=1:Table.Properties.CustomProperties.NumVars;
    else
        Table.Properties.CustomProperties.output=output;
    end

catch
    disp('Your Matlab version does not support custom properties on tables!')
    disp('Creating auxiliar table (ATable) as a file, please do not erase ATable')
    NumVars=length(out_names);
    Table2=table(NumVars);
    Table2.Domain=domain;
    Table2.rMethod=rMethod;
    Table2.Solver=solver;
    Table2.Kind=kind;
    Table2.Vars=out_names;
    Table2.Fixed=[];
     
     if isempty(output)
        Table2.output=1;
    else
        Table2.output=output;
     end
     save('ATable','Table2');
end
Table(fixed,:)=[];

%Function to manage fixed variables
    function input_control=fixing(pars,input_control,fixed)
        input_control(~fixed)=pars;
    end

end