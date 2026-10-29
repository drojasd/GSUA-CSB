classdef gsuaCoreFunctionsTest < matlab.unittest.TestCase
    %GSUACOREFUNCTIONSTEST Smoke/behavior tests for core functions that lacked direct coverage.
    %
    %   gsua_costcutoff, gsua_covmetric, gsua_sa and gsua_oatr were exercised only indirectly.
    %   These pin their core contracts so a refactor cannot quietly break them.

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

    methods (TestMethodTeardown)
        function closeFigs(~)
            close all
        end
    end

    methods (Access = private)
        function T = pkTable(~)
            T = gsua_dataprep('pkAbsorptionModel', [0.6 3.0; 0.05 0.5; 5 40], ...
                'domain', [0 24], 'names', {'ka','ke','V'}, 'out_names', {'c'});
            T.Nominal = [1.2; 0.25; 15];
        end
    end

    methods (Test)
        % ---- gsua_costcutoff ----
        function costcutoffRtolKeepsClose(testCase)
            cost = [1 1.05 1.1 10 20];
            keep = gsua_costcutoff(cost, 'rtol', 0.1, 1e-8, 3, 1);
            testCase.verifyEqual(keep, logical([1 1 1 0 0]));
        end

        function costcutoffGapCutsAtTheJump(testCase)
            cost = [1 1.02 1.05 1.1 100 200];
            keep = gsua_costcutoff(cost, 'gap', 0.1, 1e-8, 3, 1);
            testCase.verifyEqual(keep(1:4), true(1,4));
            testCase.verifyFalse(any(keep(5:6)));
        end

        function costcutoffMinKeepFloorHolds(testCase)
            % Even an aggressive cutoff keeps at least minKeep runs.
            cost = [1 100 200 300];
            keep = gsua_costcutoff(cost, 'rtol', 0.01, 1e-8, 3, 3);
            testCase.verifyGreaterThanOrEqual(sum(keep), 3);
        end

        % ---- gsua_covmetric ----
        function covmetricReturnsFiniteScoresAndOrderedBand(testCase)
            len = 15; reps = 30;
            t = linspace(0, pi, len);
            ydata = sin(t);
            rng(0); Y = ydata + 0.1*randn(reps, len);
            [cost_data, cost_band, P5, P50, P95] = gsua_covmetric(Y, ydata);
            testCase.verifyTrue(isscalar(cost_data) && isfinite(cost_data));
            testCase.verifyTrue(isscalar(cost_band) && isfinite(cost_band));
            testCase.verifyTrue(all(P5 <= P50) && all(P50 <= P95));
            testCase.verifySize(P50, [1 len]);
        end

        % ---- gsua_sa ----
        function saComputesIndicesAndShapedOutputs(testCase)
            T = testCase.pkTable();
            M = gsua_dmatrix(T, 200);
            [T2, Y, J] = gsua_sa(M, T, 'SensMethod', 'Xiao', 'parallel', false, 'verbose', false);
            testCase.verifyEqual(height(T2), 3);
            testCase.verifyTrue(all(isfinite(T2.Si)) && all(isfinite(T2.STi)));
            % Total-order mass is at least first-order mass (holds up to sampling noise).
            testCase.verifyGreaterThanOrEqual(sum(T2.STi) + 1e-6, sum(T2.Si));
            % Y is the full output matrix (one row per run), J its scalar MSE per run.
            testCase.verifySize(J, [size(Y,1) 1]);
        end

        % ---- gsua_oatr ----
        function oatrCapsUnboundedExpansion(testCase)
            % ka is insensitive on its upper side for this model; oatr used to run its
            % bound out to ~realmax (1e308). It now caps at a finite ceiling and warns.
            T = testCase.pkTable();
            Toat = testCase.verifyWarning(@() gsua_oatr(T, 'lim', 0.2), ...
                'GSUA:gsua_oatr:ExpansionCeiling');
            testCase.verifyEqual(height(Toat), 3);
            testCase.verifyTrue(all(isfinite(Toat.Range(:))), 'Ranges must be finite.');
            testCase.verifyLessThan(max(abs(Toat.Range(:))), realmax/2);
            testCase.verifyTrue(all(Toat.Range(:,1) <= Toat.Range(:,2)), ...
                'Every expanded range must have lower <= upper.');
        end
    end
end
