function [full, names] = gsua_fullpars(T, freeValues)
%GSUA_FULLPARS Reinsert fixed parameters into a free-parameter vector.
%
%   [full, names] = gsua_fullpars(T, freeValues)
%
%   A user-defined model pins a factor by giving it a degenerate range ([v v]); GSUA_DATAPREP then
%   drops that row from T entirely, so T (and every estimate produced from it, e.g. T.Estfmincon)
%   carries only the FREE factors. To evaluate a derived quantity, or to compare against the full
%   ground truth, you need the complete ordered vector with the fixed values put back. This utility
%   does exactly that: it places freeValues at the free positions and the pinned values at the
%   fixed positions, in the original factor order, and returns the matching names.
%
%   Inputs:
%     T          <-- summary table from GSUA_DATAPREP for a USER-DEFINED model (Kind 5/6). The fixed
%                    bookkeeping is read from its CustomProperties (FullNames/FixedMask/FixedValues),
%                    with a fallback to the stored `creator` for tables built before v1.18.
%     freeValues <-- the free-factor values, as a vector (row or column) of length equal to the
%                    number of free factors (height(T)), or an Np_free x K matrix whose COLUMNS are
%                    each a free vector (e.g. a whole estimate ensemble T.Estfmincon). Each column is
%                    expanded independently.
%
%   Outputs:
%     full  <-- full ordered parameter vector (Ntot x 1), or Ntot x K if freeValues had K columns.
%     names <-- Ntot x 1 cellstr of the full ordered factor names.
%
%   Example:
%     [T,~] = gsua_dataprep('myModel',[0 10; 5 5; 0 1]);   % factor 2 pinned at 5 -> dropped from T
%     [T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon');
%     [theta, names] = gsua_fullpars(T, T.Estfmincon(:,1));  % best fit, with the fixed 5 restored
%
%   See also GSUA_DATAPREP, GSUA_USERDEFINED, GSUA_DEVAL.

cp = [];
try
    cp = T.Properties.CustomProperties;
catch
end

names = [];
mask = [];
fixedVals = [];

% Preferred source: the bookkeeping recorded by gsua_userdefined at build time (v1.18+).
if ~isempty(cp)
    try
        if ~isempty(cp.FullNames)
            names     = reshape(cellstr(cp.FullNames),[],1);
            mask      = logical(reshape(cp.FixedMask,[],1));
            fixedVals = reshape(cp.FixedValues,[],1);
        end
    catch
        names = [];
    end
end

% Fallback for older tables: reconstruct from the stored creator (user-defined only).
if isempty(names)
    if isempty(cp)
        error('GSUA:gsua_fullpars:NoMetadata', ...
            'T has no CustomProperties; gsua_fullpars needs a user-defined-model table from gsua_dataprep.');
    end
    try
        cr = cp.creator;
    catch
        cr = [];
    end
    if isempty(cr) || ~isfield(cr,'names') || isempty(cr.names) || ~isfield(cr,'Range')
        error('GSUA:gsua_fullpars:Unsupported', ...
            ['gsua_fullpars supports user-defined models (built by gsua_userdefined). This table ' ...
             'lacks the FullNames/creator metadata (symbolic/Simulink models are not supported).']);
    end
    names = reshape(cellstr(cr.names),[],1);
    freeNames = reshape(T.Properties.RowNames,[],1);
    mask = ~ismember(names, freeNames);
    % creator.Range is stored before the rMethod transform; for 'range'/'normal' it equals the
    % table range, so the pinned value is Range(fixed,1). Other rMethods are not reconstructible here.
    rm = 'range';
    try, rm = cp.rMethod; catch, end
    if ~any(strcmp(rm,{'range','normal'}))
        warning('GSUA:gsua_fullpars:rMethod', ...
            'Fixed values are recovered from the raw creator range; with rMethod=''%s'' they may be off.', rm);
    end
    fixedVals = nan(numel(names),1);
    fixedVals(mask) = cr.Range(mask,1);
end

nTot  = numel(names);
nFree = nTot - sum(mask);

% Accept a vector (one free set) or a matrix whose columns are each a free set.
if isvector(freeValues)
    freeValues = reshape(freeValues,[],1);
end
if size(freeValues,1) ~= nFree
    error('GSUA:gsua_fullpars:BadLength', ...
        'freeValues must have %d rows (the number of free factors); got %d.', nFree, size(freeValues,1));
end

K = size(freeValues,2);
full = zeros(nTot,K);
full(mask,:)  = repmat(fixedVals(mask),1,K);
full(~mask,:) = freeValues;
end
