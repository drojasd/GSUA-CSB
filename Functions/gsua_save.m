function gsua_save(T,filename)
%GSUA_SAVE Save a GSUA table to disk in a portable form
%
%   gsua_save(T)
%   gsua_save(T,filename)
%
%   Removes the CustomProperties.Solver function handle from T before
%   saving, since it can reference local functions, symbolic engine
%   state, or Simulink model handles that do not serialize reliably
%   across sessions or MATLAB versions. Use gsua_load after loading the
%   file back in to reconstruct the Solver handle.
%
%   Inputs:
%     T        <-- GSUA table (from gsua_dataprep or a downstream step)
%     filename <-- (optional) output .mat file name. Default: 'portable.mat'
%
%   Example:
%     gsua_save(T,'myModel.mat')
%
%   See also GSUA_LOAD, GSUA_DATAPREP.
    if nargin<2
        filename='portable.mat';
    end
    T=rmprop(T,'Solver');
    save(filename,'T');
end