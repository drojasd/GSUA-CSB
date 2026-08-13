classdef gsuaNoisefloorWarningTest < matlab.unittest.TestCase
    %GSUANOISEFLOORWARNINGTEST Required "best fit exceeds noise floor" warning.
    %
    %   If the best pool member's cost is worse than the calibrated threshold, gsua_noisefloor
    %   must warn loudly (not pass silently) -- the model is being rejected as a description of
    %   the data. res is supplied directly (deliberately huge), which is sufficient to force this
    %   regardless of how well the underlying fit actually converged.

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
        function bestCostAboveThresholdWarns(testCase)
            [T, xdata, ydata] = noisefloorTestFixture();
            Est = T.Estfmincon;
            res = 1e12;   % arbitrarily far worse than any noise-floor threshold on this fixture

            testCase.verifyWarning( ...
                @() gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 50, 'seed', 4), ...
                'gsua_noisefloor:BestFitExceedsNoiseFloor');
        end

        function bestCostBelowThresholdDoesNotWarn(testCase)
            % The fixture's noiseless ydata also triggers the (separately-tested) quasi-Poisson/
            % global-NB degenerate-dispersion warnings regardless of res -- suppress those two
            % here so this test isolates just the one warning under test.
            warning('off', 'gsua_noisefloor:QuasiPoissonDegenerate');
            warning('off', 'gsua_noisefloor:UnderdispersedNB');
            testCase.addTeardown(@() warning('on', 'gsua_noisefloor:QuasiPoissonDegenerate'));
            testCase.addTeardown(@() warning('on', 'gsua_noisefloor:UnderdispersedNB'));

            [T, xdata, ydata] = noisefloorTestFixture();
            Est = T.Estfmincon;
            res = 1e-9;   % essentially the noiseless fit's own cost -- should clear the floor

            testCase.verifyWarningFree( ...
                @() gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 50, 'seed', 4));
        end
    end
end
