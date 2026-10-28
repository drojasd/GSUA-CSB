classdef gsuaPeSaveDefaultTest < matlab.unittest.TestCase
    %GSUAPESAVEDEFAULTTEST gsua_pe no longer writes Estimations.mat unless asked to.
    %
    %   'save' defaulted to true and wrote Estimations.mat into the current folder after every
    %   estimation, silently overwriting any file of that name. Since v1.15 it is opt-in, like
    %   the Python port's save_path. 'save' and 'timer' used to share one default variable, so
    %   these tests also pin that turning 'save' off did not turn the progress readout off.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot, 'Functions'), IncludingSubfolders=false));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot, 'Examples'), IncludingSubfolders=false));
        end
    end

    methods (TestMethodSetup)
        function runInATemporaryFolder(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
        end
    end

    methods (Access = private)
        function [T, xdata, ydata] = problem(~)
            T = gsua_dataprep('pkAbsorptionModel', [0.6 3.0; 0.05 0.5; 5 40], ...
                'domain', [0 24], 'out_names', {'c'});
            truth = [1.2; 0.25; 15]; T.Nominal = truth;
            xdata = [0.5 1 2 4 8 12 24];
            ydata = gsua_eval(truth, T, xdata, [], false, false, false);
        end
    end

    methods (Test)
        function defaultWritesNothing(testCase)
            [T, xdata, ydata] = testCase.problem();
            evalc('gsua_pe(T, xdata, ydata, ''solver'', ''lsqc'', ''N'', 1);');
            testCase.verifyFalse(isfile('Estimations.mat'));
        end

        function saveTrueStillWritesTheCheckpoint(testCase)
            [T, xdata, ydata] = testCase.problem();
            evalc('gsua_pe(T, xdata, ydata, ''solver'', ''lsqc'', ''N'', 1, ''save'', true);');
            testCase.verifyTrue(isfile('Estimations.mat'));
        end

        function timerStaysOnByDefault(testCase)
            % Guards the decoupling: the per-estimation readout must still print by default.
            [T, xdata, ydata] = testCase.problem();
            out = evalc('gsua_pe(T, xdata, ydata, ''solver'', ''lsqc'', ''N'', 1);');
            testCase.verifySubstring(out, 'Estimation 1');
        end
    end
end
