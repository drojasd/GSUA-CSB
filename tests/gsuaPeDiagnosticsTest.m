classdef gsuaPeDiagnosticsTest < matlab.unittest.TestCase
    %GSUAPEDIAGNOSTICSTEST gsua_pe fails loudly, not silently, when estimation cannot work.
    %
    %   A model that does not evaluate to finite output used to leave gsua_pe returning an
    %   all-Inf res with no cause (the reported "fcosto = 1x50 todo Inf"). gsua_pe now runs a
    %   pre-flight evaluation at the nominal and errors with GSUA:gsua_pe:ModelNotFinite, and a
    %   post-loop guard errors with GSUA:gsua_pe:AllEstimationsFailed if every estimation still
    %   came back non-finite. A healthy model must of course still estimate without error.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot, 'Functions'), IncludingSubfolders=false));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot, 'Examples'), IncludingSubfolders=false));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                testDir, IncludingSubfolders=false));
        end
    end

    methods (Test)
        function nonFiniteModelErrorsEarly(testCase)
            % The pre-flight check must fire before any optimizer runs.
            T = gsua_dataprep('peNonFiniteModel', [0 1], 'domain', [0 10], 'out_names', {'y'});
            xdata = linspace(0, 10, 20);
            ydata = ones(1, 20);
            testCase.verifyError(@() gsua_pe(T, xdata, ydata, 'solver', 'fmincon', ...
                'N', 2, 'save', false, 'timer', false), 'GSUA:gsua_pe:ModelNotFinite');
        end

        function healthyModelStillEstimates(testCase)
            % A finite model must not trip the new guards.
            ranges = [0.6 3.0; 0.05 0.5; 5 40];
            T = gsua_dataprep('pkAbsorptionModel', ranges, 'domain', [0 24], ...
                'names', {'ka','ke','V'}, 'out_names', {'c'});
            truth = [1.2; 0.25; 15]; T.Nominal = truth;
            xdata = [0.5 1 2 4 8 12 24];
            ydata = gsua_eval(truth, T, xdata, [], false, false, false);
            [~, res] = gsua_pe(T, xdata, ydata, 'solver', 'lsqc', 'N', 2, ...
                'save', false, 'timer', false);
            testCase.verifyTrue(all(isfinite(res)));
        end
    end
end
