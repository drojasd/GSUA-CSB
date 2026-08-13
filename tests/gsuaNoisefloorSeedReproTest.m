classdef gsuaNoisefloorSeedReproTest < matlab.unittest.TestCase
    %GSUANOISEFLOORSEEDREPROTEST Reproducibility via 'seed', and caller rng state restoration.

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
        function sameSeedGivesBitIdenticalResults(testCase)
            [T, xdata, ydata] = noisefloorTestFixture();
            Est = T.Estfmincon;
            res = 0.05;

            out1 = gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 100, 'seed', 7);
            out2 = gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 100, 'seed', 7);

            testCase.verifyEqual(out1.cost, out2.cost);
            testCase.verifyEqual(out1.threshold, out2.threshold);
            testCase.verifyEqual(out1.accepted, out2.accepted);
        end

        function differentSeedsGiveDifferentDraws(testCase)
            [T, xdata, ydata] = noisefloorTestFixture();
            Est = T.Estfmincon;
            res = 0.05;

            out1 = gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 100, 'seed', 7);
            out2 = gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 100, 'seed', 8);

            testCase.verifyNotEqual(out1.cost, out2.cost);
        end

        function callerGlobalRngStateIsRestored(testCase)
            [T, xdata, ydata] = noisefloorTestFixture();
            Est = T.Estfmincon;
            res = 0.05;

            rng(42, 'twister');
            stateBefore = rng;
            gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 50, 'seed', 9);
            stateAfter = rng;

            testCase.verifyEqual(stateAfter.Seed, stateBefore.Seed);
            testCase.verifyEqual(stateAfter.Type, stateBefore.Type);
            testCase.verifyEqual(stateAfter.State, stateBefore.State);
        end
    end
end
