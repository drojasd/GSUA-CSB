classdef gsuaPeSolverSmokeTest < matlab.unittest.TestCase
    %GSUAPESOLVERSMOKETEST Every gsua_pe solver runs without error.
    %
    %   Guards the solver-dispatch calls in gsua_pe -- in particular the previously-broken
    %   'psearch' (passed Np instead of a start point) and 'annealing' (missing x0, shifting every
    %   argument). Each solver is exercised once on a tiny SIR problem; a solver whose optimizer is
    %   not installed is skipped (assumption), not failed.

    properties (TestParameter)
        solver = {'lsqc','lsqn','ga','particle','psearch','surrogate','annealing','fmincon'};
    end

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
        function solverRunsWithoutError(testCase,solver)
            fnMap = struct('lsqc','lsqcurvefit','lsqn','lsqnonlin','ga','ga', ...
                'particle','particleswarm','psearch','patternsearch','surrogate','surrogateopt', ...
                'annealing','simulannealbnd','fmincon','fmincon');
            reqFn = fnMap.(solver);
            testCase.assumeNotEmpty(which(reqFn), ...
                sprintf('%s not available on this installation; skipping %s.', reqFn, solver));

            ranges = [1e3 1e4; 0 10; 10 100; 0 1; 0 1; 0 1];
            [T,~] = gsua_dataprep('sir_sencillo', ranges, 'domain',[1 200], ...
                'names',{'S0','I0','R0','beta','gamma','alpha'}, 'out_names',{'S','I','R'});
            T.Properties.CustomProperties.output = 2;
            truth = [5000;5;50;0.5;0.2;0.01]; xs = 1:8;
            ys = squeeze(gsua_eval(truth,T,xs,[],false,false,false)); ys = ys(:).';

            rng(0);
            warnState = warning('off','all');
            cleanup = onCleanup(@() warning(warnState)); %#ok<NASGU>
            [~,r] = gsua_pe(T,xs,ys,'solver',solver,'N',1,'margin',0.3,'save',false,'timer',false);
            testCase.verifyNumElements(r,1);
            testCase.verifyTrue(isfinite(r(1)), sprintf('%s returned non-finite res.', solver));
        end
    end
end
