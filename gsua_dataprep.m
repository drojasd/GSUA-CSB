function [T,out1]=gsua_dataprep(p1,p2,p3,p4,varargin)
% Consolidation of data and environment for (c) GSUA-RC Toolbox
%
% Working with simulink models
% [T,tout]=gsua_dataprep(model_name,Ranges,Parameter_names)
% Parameters:
% model_name      <-- file name of the simulink model
% Ranges          <-- double array of Npx2
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
% Working with User-defined functions, matlab
% [T,solver]=gsua_dataprep(Function,Ranges)
% Parameters:
% Function       <-- Is the name of the user-defined function (char)
% Ranges       <-- double array of Npx2
% Outputs:
% T      <-- table array with all required system information
% solver <-- handle function to solve the system
%
% To explore additional configuration, open the user guide:
% open gsua_userguide <--- User-defined functions userguide

if isstring(p1)||ischar(p1)
    type=exist(p1,'file');
    switch type
        case 4
            disp('Setting environment to work with simulink')
            if nargin<4
                [T,out1] = sens_dataprep(p1,p2,p3,varargin{:});
            else
                [T,out1] = sens_dataprep(p1,p2,p3,p4,varargin{:});
            end        
        case 2
            disp('Setting environment to work with user-defined function')
            if nargin==2
                [T,out1] = gsua_userdefined(p1,p2);
            else
                [T,out1] = gsua_userdefined(p1,p2,p3,p4,varargin{:});
            end
        otherwise
            disp(p1,' Is not recognized, check if it is on your current path')
            return
        
    end

else
    disp('Setting environment to work with symbolic Matlab')
    [T,out1] = gsua_dpmat(p1,p2,p3,p4,varargin{:});
end

end