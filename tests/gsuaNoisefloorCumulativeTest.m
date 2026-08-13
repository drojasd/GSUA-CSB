classdef gsuaNoisefloorCumulativeTest < matlab.unittest.TestCase
    %GSUANOISEFLOORCUMULATIVETEST Cumulative-output handling and real-data NaN mask replication.
    %
    %   Row 2 of the fixture model is treated as cumulative. A deliberate NaN gap is punched into
    %   its real ydata; the synthetic bootstrap sample's NaN pattern for that row must match
    %   ydata's own exactly (verified via out.exampleSynthetic, since gsua_noisefloor otherwise
    %   only returns aggregate cost distributions). No NaN may leak into the final cost values.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot, 'Functions'), IncludingSubfolders=false));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                testDir, IncludingSubfolders=false));
        end
    end

    methods (Test)
        function syntheticNanMaskMatchesRealDataInCumulativeSpace(testCase)
            [T, xdata, ydata] = noisefloorTestFixture();
            ydata(2, 6:8) = NaN;   % gap in the cumulative row

            Est = T.Estfmincon;
            res = 0.01;   % small, non-degenerate cost so dispersion estimation stays well-posed

            out = gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 50, 'seed', 2, ...
                'cumulative', [false true]);

            synthRow2 = out.exampleSynthetic.quasipoisson(2, :);
            testCase.verifyTrue(all(isnan(synthRow2(6:8))));
            testCase.verifyTrue(all(isfinite(synthRow2([1:5, 9:end]))));

            % Row 1 (non-cumulative, no real gap) must be fully finite throughout.
            synthRow1 = out.exampleSynthetic.quasipoisson(1, :);
            testCase.verifyTrue(all(isfinite(synthRow1)));

            testCase.verifyTrue(all(isfinite(out.cost)));
            testCase.verifyFalse(any(isnan(out.cost)));
        end

        function cumulativeRowSyntheticSampleIsNondecreasing(testCase)
            % A row treated as cumulative is reconstructed via cumsum of non-negative incident
            % draws -- its synthetic sample must be non-decreasing outside the NaN gap, confirming
            % the cumulative re-accumulation actually happened rather than treating the row as
            % flat incident data.
            [T, xdata, ydata] = noisefloorTestFixture();

            out = gsua_noisefloor(T, xdata, ydata, T.Estfmincon, 0.01, 'B', 50, 'seed', 3, ...
                'cumulative', [false true]);

            synthRow2 = out.exampleSynthetic.quasipoisson(2, :);
            testCase.verifyTrue(all(diff(synthRow2) >= 0));
        end
    end
end
