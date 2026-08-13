classdef gsuaNoisefloorQuasiPoissonFallbackTest < matlab.unittest.TestCase
    %GSUANOISEFLOORQUASIPOISSONFALLBACKTEST Quasi-Poisson degenerate-dispersion fallback.
    %
    %   With noiseless data (residuals identically zero at theta*), the pooled Pearson dispersion
    %   phi is <= 1 by construction (num=0), which is degenerate for the quasi-Poisson
    %   reparameterization (k = mu/(phi-1) requires phi>1) -- gsua_noisefloor must fall back to
    %   plain Poisson draws with a warning, not error or silently produce a nonsensical k.

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
        function noiselessFitTriggersQuasiPoissonFallback(testCase)
            [T, xdata, ydata] = noisefloorTestFixture();
            Est = T.Estfmincon;
            res = 0;   % noiseless fit converges to (near) exactly zero cost

            testCase.verifyWarning( ...
                @() gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 2000, 'seed', 1), ...
                'gsua_noisefloor:QuasiPoissonDegenerate');

            warning('off', 'gsua_noisefloor:QuasiPoissonDegenerate');
            warning('off', 'gsua_noisefloor:UnderdispersedNB');
            cleanupObj = onCleanup(@() warning('on', 'all'));
            out = gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 2000, 'seed', 1);

            testCase.verifyTrue(out.quasiFallback);
            % The fallback draws from a different point in the shared RNG stream than the
            % separately-run 'poisson' model entry (each model's bootstrap runs in sequence off
            % one seeded stream, not independently re-seeded), so the two B x 1 distributions are
            % NOT bit-identical -- but both are literally Poisson(mu) draws, so their summary
            % statistics at B=2000 should agree closely.
            testCase.verifyEqual(mean(out.byModel.quasipoisson.cost), mean(out.byModel.poisson.cost), ...
                'RelTol', 0.15);
            testCase.verifyEqual(std(out.byModel.quasipoisson.cost), std(out.byModel.poisson.cost), ...
                'RelTol', 0.3);
        end
    end
end
