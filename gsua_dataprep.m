function [T,out1]=gsua_dataprep(p1,p2,p3,p4,varargin)
% Consolidation of data and environment for (c) GSUA-RC Toolbox
%
% Working with simulink models
% [T,tout]=gsua_dataprep(model_name,Ranges,Parameter_names)
% Parameters:
% model_name      <-- file name of the simulink model
% Ranges          <-- cell array of Npx2
% Parameter_names <-- cell array with the parameter names of the model
% Outputs:
% T    <-- table array with all required system information
% tout <-- array of output time for the model
%
% Working with symbolic toolbox, matlab
% [T,solver]=gsua_dataprep(odes,vars,domain,model_name)
% Parameters:
% odes       <-- array with the system of differential equations
% vars       <-- array of state variables for the system of equations
% domain     <-- array with the domain of the system of equations
% model_name <-- name for the output file of the system
% Outputs:
% T      <-- table array with all required system information
% solver <-- handle function for ode solvers
%
% To explore additional configuration, open the user guide:
% open gsua_userguide <--- Matlab userguide

if isstring(p1)||ischar(p1)
    if nargin<4
        [T,out1] = sens_dataprep(p1,p2,p3,varargin{:});
    else
        [T,out1] = sens_dataprep(p1,p2,p3,p4,varargin{:});
    end
    disp('Setting environment to work with simulink')
    

else
    disp('Setting environment to work with symbolic Matlab')
    [T,out1] = gsua_dpmat(p1,p2,p3,p4,varargin{:});
end

end