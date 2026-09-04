classdef gsuaDmatrixJointSeedTest < matlab.unittest.TestCase
    %GSUADMATRIXJOINTSEEDTEST Seeded joint sampling is reproducible and RNG-hygienic.
    %
    %   Bands built for a paper have to be reproducible, so 'Seed' must fully determine the
    %   draw. It must also not leak: seeding this function is a local request, and silently
    %   moving the caller's global RNG stream would make every later random call in the
    %   caller's script depend on whether a band happened to be sampled.

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
        function sameSeedReproducesExactly(testCase)
            T = jointTestFixture();
            for jt = {'Bootstrap','SmoothBootstrap','Gaussian'}
                M1 = gsua_dmatrix(T, 200, 'Method','Joint', 'JointType', jt{1}, 'Seed', 2024);
                M2 = gsua_dmatrix(T, 200, 'Method','Joint', 'JointType', jt{1}, 'Seed', 2024);
                testCase.verifyEqual(M1, M2);
            end
        end

        function differentSeedsGiveDifferentDraws(testCase)
            T = jointTestFixture();
            M1 = gsua_dmatrix(T, 200, 'Method','Joint', 'Seed', 1);
            M2 = gsua_dmatrix(T, 200, 'Method','Joint', 'Seed', 2);
            testCase.verifyNotEqual(M1, M2);
        end

        function callerGlobalRngStateIsRestored(testCase)
            T = jointTestFixture();
            rng(42, 'twister');
            before = rng;
            gsua_dmatrix(T, 100, 'Method','Joint', 'Seed', 7);
            after = rng;
            testCase.verifyEqual(after.Seed, before.Seed);
            testCase.verifyEqual(after.Type, before.Type);
            testCase.verifyEqual(after.State, before.State);
        end

        function unseededCallsAdvanceTheStreamNormally(testCase)
            % Without 'Seed' the function must behave like every other sampler in the
            % toolbox: draw from the caller's stream and leave it advanced.
            T = jointTestFixture();
            rng(11, 'twister');
            M1 = gsua_dmatrix(T, 100, 'Method','Joint');
            M2 = gsua_dmatrix(T, 100, 'Method','Joint');
            testCase.verifyNotEqual(M1, M2);
        end
    end
end
