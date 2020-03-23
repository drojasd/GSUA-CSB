function [ParT,tOut] = sens_dataprep(model,Ranges,ParIn,varargin)
% Consolidation of data and simulink environment for (c) GSUA-RC Toolbox
%
% [T,tout]=sens_dataprep(model_name,Ranges,Parameter_names)
% Parameters:
% model_name      <-- file name of the simulink model
% Ranges          <-- cell array of Npx2
% Parameter_names <-- cell array with the parameter names of the model
% Outputs:
% T    <-- table array with all required system information
% tout <-- array of output time for the model
%
% To explore additional configuration, open the user guide:
% open gsua_simuserguide <--- simulink userguide

%% Definition of default values for SampleMethod,rMethod,nominals,additional

p=inputParser;

defaultrMethod='range';
validrMethod={'percent','range','normal', 'std'};
checkrMethod = @(x) any(validatestring(x,validrMethod));

defaultTypeModel='Dynamic';
validTypeModel={'Dynamic','Static'};
checkTypeModel = @(x) any(validatestring(x,validTypeModel));

defaultNominals={};
defaultStep=[];

checkPar= @(x) any([istable(x),iscell(x)]);
checkRanges = @(x) any([isnumeric(x),iscell(x)]);

addRequired(p,'model',@ischar);
addRequired(p,'Ranges',checkRanges);
addRequired(p,'ParIn',checkPar);
addParameter(p,'rMethod',defaultrMethod,checkrMethod);
addParameter(p,'modelkind',defaultTypeModel,checkTypeModel);
addParameter(p,'nominal',defaultNominals,checkRanges);
addParameter(p,'Step',defaultStep,@isnumeric);

parse(p,model,Ranges,ParIn,varargin{:})

Ranges=p.Results.Ranges;

Np = size(Ranges,1);
if isa(Ranges,'double')
    Ranges=mat2cell(Ranges,repelem(1,Np),[1 1]);
end
ParOut=cell(Np,4);
  
if nargout~= 1
%% Creating temporal ParOut
model=p.Results.model;
ParIn=p.Results.ParIn;
rMethod=p.Results.rMethod;
TypeModel=p.Results.modelkind;

ParOut(:,1)=ParIn;

if size(Ranges,2)>2
    ParOut(:,2)=Ranges(:,3);
else
    ParOut(:,2)=repelem(cellstr(rMethod),Np)';
end


%% Ranges calculation

for k = 1:Np
    switch ParOut{k,2}
        case 'std'%std
            ParOut{k,4} = Ranges{k,1};
            ParOut{k,3} = [Ranges{k,1}-Ranges{k,2} ; Ranges{k,1}+Ranges{k,2}];
        case 'percent'%percent
            ParOut{k,4} = Ranges{k,1};
            ParOut{k,3} = [Ranges{k,1}-Ranges{k,1}*Ranges{k,2}/100 ; Ranges{k,1}+Ranges{k,1}*Ranges{k,2}/100];
        case 'range'%range
            ParOut{k,4} = (Ranges{k,1}+Ranges{k,2})/2;
            ParOut{k,3} = [Ranges{k,1} ; Ranges{k,2}];
        case 'normal'
            ParOut{k,4} = Ranges{k,1};
            ParOut{k,3} = Ranges{k,2};
    end
end

if isempty(p.Results.nominal)
    disp('Setting nominal values')
else
    if isa(p.Results.nominal,'double')
        ParOut(:,4)=mat2cell(p.Results.nominal',repelem(1,Np),1);
    else
        ParOut(:,4)=p.Results.nominal;
    end
end

%% Model configuration
load_system(model)
set_param(model,'SaveFormat','Array'); %Format configuration
set_param(model,'SaveOutput','on','SaveTime','on'); %Outputs warranted
set_param(model,'OutputSaveName','yout','OutputTimes','tout','ReturnWorkspaceOutputs', 'on'); %Output names configuration

if ~isempty(p.Results.Step) %If the user specifies a step, use it
    set_param(model,'Fixedstep',mat2str(p.Results.Step));
end

switch TypeModel 
    case 'Static' %For static models the initial and final time are the same
        set_param(model,'StartTime','0');
        set_param(model, 'StopTime','0')
    case 'Dynamic'
    otherwise
        disp('Unknown model type')
        return
end

hws = get_param(model, 'ModelWorkspace');

for i=1:Np
    hws.assignin(ParOut{i,1}, ParOut{i,4});
end

simOutArray = sim(model);
%y_Nom = simOutArray.get('yout')';
tOut = simOutArray.get('tout');

switch TypeModel
    case'Dynamic' 
        if isempty(p.Results.Step)
            set_t=mat2str(tOut(2)-tOut(1));
            set_param(model,'Fixedstep',set_t);
        end
end

close_system(model, 1);
else
    ParOut(:,3)=mat2cell([Ranges{:,1};Ranges{:,2}]',repelem(1,Np),2);
end

if nargout ~=1
ParT=table([ParOut{:,3}]',[ParOut{:,4}]','RowNames',ParOut(:,1),'VariableNames',{'Range','Nominal'});
end
try
 ParT=addprop(ParT,{'Kind','rMethod','Fixed','output','Vars','Domain','Tout'},{'table','table','table','table','table','table','table'});
 ParT.Properties.CustomProperties.Kind=1;
 ParT.Properties.CustomProperties.rMethod=rMethod;
 ParT.Properties.CustomProperties.Fixed=[];
 ParT.Properties.CustomProperties.output=1;
 ParT.Properties.CustomProperties.Vars={'Output'};
 ParT.Properties.CustomProperties.Domain=[min(tOut),max(tOut)];
 ParT.Properties.CustomProperties.Tout=tOut;
disp('All done!')
catch
    disp('Your Matlab version does not support custom properties on tables!')
    disp('Creating auxiliar table (ATable) as a file, please do not erase ATable')
    Table2=table(model);
    Table2.Kind=1;
    Table2.Fixed=[];
    Table2.output=1;
    Table2.Vars={'Output'};
    Table2.Domain=[min(tOut),max(tOut)];
    Table2.Tout=tOut;
    save('ATable','Table2');
    disp('All done!')
end
end
