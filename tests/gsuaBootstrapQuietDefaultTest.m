classdef gsuaBootstrapQuietDefaultTest < matlab.unittest.TestCase
    %GSUABOOTSTRAPQUIETDEFAULTTEST gsua_bootstrap runs its refits quietly by default (v1.17.1 perf).
    %
    %   With no 'opt', the B internal refits used to call fmincon with default options, which prints a
    %   "Local minimum found" block per start; on an interactive desktop that console rendering turned
    %   a ~54 s bootstrap into ~234 s. gsua_bootstrap now builds Display-off options for the solver by
    %   default. Display is cosmetic, so the ensemble for a given seed must be unchanged -- these
    %   tests check both: (1) the default path produces no optimizer chatter, and (2) it is identical
    %   to passing explicit Display-off options.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot,'Functions'), IncludingSubfolders=false));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot,'Examples'), IncludingSubfolders=false));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                testDir, IncludingSubfolders=false));
        end
    end

    methods (Test)
        function defaultIsQuietAndResultUnchanged(testCase)
            [T,xdata,ydata] = sirSencilloFixture();
            opt = optimoptions('fmincon','Display','off');
            rng(2);
            [T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.3,'N',6, ...
                'opt',opt,'save',false,'timer',false);

            % (1) The default (no 'opt') path must not spill optimizer output.
            chatter = evalc(['Edef = gsua_bootstrap(T,xdata,ydata,T.Estfmincon,res,' ...
                '''B'',3,''solver'',''fmincon'',''N'',3,''seed'',5,''timer'',false);']);
            testCase.verifyEmpty(regexp(chatter,'Local minimum','once'), ...
                'The default refits must run with Display off (no fmincon chatter).');

            % (2) Same seed, explicit Display-off options -> identical ensemble (Display is cosmetic).
            Eexp = gsua_bootstrap(T,xdata,ydata,T.Estfmincon,res,'B',3,'solver','fmincon', ...
                'N',3,'opt',opt,'seed',5,'timer',false);
            testCase.verifyEqual(Edef,Eexp, ...
                'The quiet default must not change the estimates for a given seed.');
        end
    end
end
