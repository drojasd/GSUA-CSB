classdef gsuaMarginAndClusterFlagTest < matlab.unittest.TestCase
    %GSUAMARGINANDCLUSTERFLAGTEST Unified margin convention (#3) and clustering-not-evaluated flag (#10).

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
        function likelihoodMarginIsRawAndWarnsOnOldConvention(testCase)
            % margin is now the RAW relative noise (0.08 = 8%); a value >= 1 looks like the pre-v1.18
            % offset convention and must warn.
            T = gsua_userdefined('noisefloorFixtureFunc',[0.1 20;0.01 2;0.01 5]);
            T.Nominal=[4;0.4;0.8];
            xs = linspace(0,4,15);
            raw  = @() gsua_likelihood(T, xs, [], 0.95, 0.1, 0.08, 0.5, 0.5, 1, 1, false, false, false, 2, false);
            old  = @() gsua_likelihood(T, xs, [], 0.95, 0.1, 1.08, 0.5, 0.5, 1, 1, false, false, false, 2, false);
            testCase.verifyWarningFree(raw, 'A raw margin (0.08) must not warn.');
            testCase.verifyWarning(old, 'GSUA:gsua_likelihood:MarginConvention');
        end

        function iaFlagsClusteringNotEvaluated(testCase)
            T = gsua_userdefined('noisefloorFixtureFunc',[0.1 20;0.01 2;0.01 5]);
            T.Nominal=[4;0.4;0.8];
            rng(0); Test = repmat(T.Nominal,1,2) + 0.01*randn(3,2);   % only 2 runs -> cannot cluster
            % cluster requested but not enough runs: warns with an id AND flags it in the output.
            f = @() gsua_ia(T, Test, false, false, false, true);
            testCase.verifyWarning(f, 'GSUA:gsua_ia:FewRunsNoCluster');
            ws = warning('off','all'); c = onCleanup(@() warning(ws));
            [~,ci] = gsua_ia(T, Test, false, false, false, true);
            testCase.verifyTrue(ci.ClusterRequested);
            testCase.verifyFalse(ci.ClusteringEvaluated);
            testCase.verifyNotEmpty(ci.Note);
        end
    end
end
