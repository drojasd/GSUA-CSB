function T=gsua_load(T)
%GSUA_LOAD Rebuild a GSUA table's Solver handle after loading from disk
%
%   T = gsua_load(T)
%
%   gsua_save strips the CustomProperties.Solver function handle before
%   saving a GSUA table (function handles to model/ODE definitions are
%   not always safe or compact to serialize). gsua_load reverses that:
%   it reads the table's stored 'creator' metadata (the original
%   arguments used to build the table) and re-runs the matching
%   gsua_dataprep/gsua_dpmat/gsua_userdefined call to reconstruct the
%   Solver handle and any other derived CustomProperties.
%
%   Inputs:
%     T <-- GSUA table previously loaded from a .mat file saved with
%           gsua_save (must still carry its CustomProperties.Kind and
%           CustomProperties.creator fields)
%
%   Outputs:
%     T <-- same table with CustomProperties.Solver (and, for Simulink-
%           or symbolic-derived tables, Domain) restored
%
%   Example:
%     load('portable.mat','T')
%     T = gsua_load(T);
%
%   See also GSUA_SAVE, GSUA_DATAPREP.

prop=T.Properties.CustomProperties;

Kind = prop.Kind;

switch Kind
    case {5,6}%userdefined
        func=prop.creator.func;
        Range=prop.creator.Range;
        try
            names=prop.creator.names;
        catch
            names=[];
        end
        nominal=prop.creator.nominal;
        vectorized=prop.creator.vectorized;
        domain=prop.Domain;
        rMethod=prop.rMethod;
        output=prop.output;
        vars=prop.Vars;
        opt=prop.copt;



        T2=gsua_dataprep(func,Range,'names',names,'domain',domain,'rMethod',rMethod,...
            'nominal',nominal,'output',output,'out_names',vars,'vectorized',vectorized,...
            'opt',opt);
        T.Properties.CustomProperties=T2.Properties.CustomProperties;
        
    case {3,4}%dpmat
        
        odes = prop.creator.odes;
        vars = prop.creator.vars;
        domain = prop.Domain;
        modelname = prop.creator.modelname;
        
        T2=gsua_dataprep(odes,vars,domain,modelname);
        T.Properties.CustomProperties.Solver=T2.Properties.CustomProperties.Solver;
        
    case 1%simulink
        
        Mname = prop.MName;
        Ranges = prop.creator.Ranges;
        ParIn = prop.creator.ParIn;
        
        T2=gsua_dataprep(Mname,Ranges,ParIn);
        T.Properties.CustomProperties.Solver=T2.Properties.CustomProperties.Solver;
        
    case 2%dpmat or userdefined
        
        try
            odes = prop.creator.odes;
            vars = prop.creator.vars;
            domain = prop.Domain;
            modelname = prop.creator.modelname;

            T2=gsua_dataprep(odes,vars,domain,modelname);
            T.Properties.CustomProperties.Solver=T2.Properties.CustomProperties.Solver;

        catch
            func=prop.creator.func;
            Range=prop.creator.Range;
            try
                names=prop.creator.names;
            catch
                names=[];
            end
            nominal=prop.creator.nominal;
            vectorized=prop.creator.vectorized;
            domain=prop.Domain;
            rMethod=prop.rMethod;
            output=prop.output;
            vars=prop.Vars;
            opt=prop.copt;



            T2=gsua_dataprep(func,Range,'names',names,'domain',domain,'rMethod',rMethod,...
                'nominal',nominal,'output',output,'out_names',vars,'vectorized',vectorized,...
                'opt',opt);
            T.Properties.CustomProperties=T2.Properties.CustomProperties;
        end
        
    otherwise
end

end