classdef gsuaPeYdataOrientTest < matlab.unittest.TestCase
    %GSUAPEYDATAORIENTTEST gsua_pe orients ydata to outputs x time for every solver.
    %
    %   Regression for the orientation bug: the gradient default FMINCON and the metaheuristics do
    %   not auto-transpose ydata the way lsqc/lsqn do, so a single-output ydata passed as a COLUMN
    %   (len x 1) used to be read as len one-point "outputs" and crash GSUA_COSTF with an opaque
    %   'badsubscript' ("Every estimation failed"). gsua_pe now orients ydata (outputs x time)
    %   up front, so a column and a row give the same result. The solver smoke test passes a row,
    %   so it never exercised this path.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot,'Functions'), IncludingSubfolders=false));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(testDir, IncludingSubfolders=false));
        end
    end

    methods (Static)
        function [T,xs,ysRow] = makeProblem()
            ranges = [1e3 1e4; 0 10; 10 100; 0 1; 0 1; 0 1];
            [T,~] = gsua_dataprep('sir_sencillo', ranges, 'domain',[1 200], ...
                'names',{'S0','I0','R0','beta','gamma','alpha'}, 'out_names',{'S','I','R'});
            T.Properties.CustomProperties.output = 2;              % single observed output
            truth = [5000;5;50;0.5;0.2;0.01]; xs = 1:8;
            ysRow = squeeze(gsua_eval(truth,T,xs,[],false,false,false));
            ysRow = ysRow(:).';                                    % 1 x nt (outputs x time)
        end
    end

    methods (Test)
        function fminconColumnYdataRuns(testCase)
            [T,xs,ysRow] = gsuaPeYdataOrientTest.makeProblem();
            w = warning('off','all'); c = onCleanup(@() warning(w)); %#ok<NASGU>
            rng(0);
            [~,r] = gsua_pe(T,xs,ysRow(:),'solver','fmincon','N',1,'margin',0.3,'save',false,'timer',false);
            testCase.verifyTrue(isfinite(r(1)), ...
                'fmincon failed with a column ydata (orientation not handled).');
        end

        function rowAndColumnAgree(testCase)
            [T,xs,ysRow] = gsuaPeYdataOrientTest.makeProblem();
            w = warning('off','all'); c = onCleanup(@() warning(w)); %#ok<NASGU>
            rng(0); [~,rRow] = gsua_pe(T,xs,ysRow,   'solver','fmincon','N',1,'margin',0.3,'save',false,'timer',false);
            rng(0); [~,rCol] = gsua_pe(T,xs,ysRow(:),'solver','fmincon','N',1,'margin',0.3,'save',false,'timer',false);
            testCase.verifyEqual(rCol(1), rRow(1), 'RelTol', 1e-9, ...
                'Row and column ydata give different costs; orientation not applied consistently.');
        end
    end
end
