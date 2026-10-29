classdef gsuaPardevalFailWarningTest < matlab.unittest.TestCase
    %GSUAPARDEVALFAILWARNINGTEST gsua_pardeval warns when model evaluations fail.
    %
    %   gsua_pardeval ran each design-matrix row inside a try/catch with an EMPTY catch, so a
    %   simulation that errored was silently left as its preallocated zeros and the batch (used
    %   by gsua_ua/gsua_sa) came back quietly incomplete. It now counts the failures and warns
    %   with GSUA:gsua_pardeval:FailedSimulations, while still returning the partial result.

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
        function failedSimulationsWarn(testCase)
            T = gsua_dataprep('pardevalErroringModel', [0 1; 0 1], 'domain', [0 10], ...
                'out_names', {'y'});
            M = gsua_dmatrix(T, 4);
            xdata = linspace(0, 10, 12);
            y = testCase.verifyWarning(@() gsua_pardeval(M, T, xdata, false), ...
                'GSUA:gsua_pardeval:FailedSimulations');
            % The batch still returns, shaped as usual, with the failed rows left as zeros.
            testCase.verifySize(y, [4, numel(xdata)]);
            testCase.verifyEqual(y, zeros(4, numel(xdata)));
        end

        function healthyBatchDoesNotWarn(testCase)
            T = gsua_dataprep('pkAbsorptionModel', [0.6 3.0; 0.05 0.5; 5 40], ...
                'domain', [0 24], 'out_names', {'c'});
            M = gsua_dmatrix(T, 4);
            xdata = [1 2 4 8 12 24];
            testCase.verifyWarningFree(@() gsua_pardeval(M, T, xdata, false));
        end
    end
end
