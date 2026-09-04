classdef gsuaDmatrixJointCorrelationTest < matlab.unittest.TestCase
    %GSUADMATRIXJOINTCORRELATIONTEST Joint sampling must carry the pool's correlation through.
    %
    %   The whole point of the joint modes is that the correlation structure among parameters
    %   survives sampling, where the marginal methods destroy it. These tests pin both halves
    %   of that claim -- and one boundary of it: the Gaussian mode reproduces the correlation
    %   COEFFICIENT well while still leaving a curved manifold (gsuaDmatrixJointManifoldTest),
    %   which is why correlation alone is not the acceptance criterion for this feature.

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
        function bootstrapPreservesPoolCorrelation(testCase)
            [T,pool] = jointTestFixture();
            M = gsua_dmatrix(T, 4000, 'Method','Joint', 'Seed', 21);
            testCase.verifyEqual(corr(M(:,1),M(:,2)), ...
                corr(pool(1,:)',pool(2,:)'), 'AbsTol', 0.05);
        end

        function marginalSamplingDestroysCorrelation(testCase)
            % The documented failure this feature exists to fix.
            T = jointTestFixture();
            M = gsua_dmatrix(T, 4000, 'Seed', 21);
            testCase.verifyLessThan(abs(corr(M(:,1),M(:,2))), 0.15);
        end

        function gaussianAlsoPreservesLinearCorrelation(testCase)
            [T,pool] = jointTestFixture();
            M = gsua_dmatrix(T, 4000, 'Method','Joint', 'JointType','Gaussian', 'Seed', 21);
            testCase.verifyEqual(corr(M(:,1),M(:,2)), ...
                corr(pool(1,:)',pool(2,:)'), 'AbsTol', 0.05);
        end

        function bootstrapInheritsPoolSpreadNotTableRange(testCase)
            % gsua_ia leaves T.Range as a CI of the MEDIAN, which is narrower than the pool
            % and narrows further as the pool grows. Joint draws must reflect the pool's own
            % spread instead -- that is the second error this mode fixes, and it only holds
            % because nothing clips the draws back to T.Range by default.
            [T,pool] = jointTestFixture();
            T.Range(1,:) = [1.4 1.6];       % a deliberately too-narrow "CI of the median"
            M = gsua_dmatrix(T, 2000, 'Method','Joint', 'Pool', pool, 'Seed', 33);
            testCase.verifyEqual(std(M(:,1)), std(pool(1,:)), 'RelTol', 0.15);
            testCase.verifyGreaterThan(max(M(:,1)), 1.6);
            testCase.verifyLessThan(min(M(:,1)), 1.4);
        end

        function clipIsAppliedWhenRequested(testCase)
            [T,pool] = jointTestFixture();
            bounds = [1.0 2.0; 0.3 3; 0.1 1.2];
            M = gsua_dmatrix(T, 500, 'Method','Joint', 'Pool', pool, ...
                'Clip', bounds, 'Seed', 33);
            testCase.verifyGreaterThanOrEqual(min(M(:,1)), 1.0);
            testCase.verifyLessThanOrEqual(max(M(:,1)), 2.0);
        end
    end
end
