classdef gsuaDmatrixJointGuardTest < matlab.unittest.TestCase
    %GSUADMATRIXJOINTGUARDTEST A pool too small to trust must be refused, not used.
    %
    %   Correlation estimated from very few points is degenerate rather than merely noisy --
    %   two points always correlate at exactly +-1 -- so joint sampling from such a pool
    %   produces a band that is confidently wrong instead of approximately right. That is a
    %   worse failure than the one this feature fixes, because it looks like success. This is
    %   not hypothetical: a real PEtab case in the project's own experiments (Boehm 2014)
    %   yields a two-point dominant cluster after clustering.
    %
    %   The guard fires on pool size and falls back to marginal sampling with a warning.

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
        function tinyPoolWarns(testCase)
            [T,pool] = jointTestFixture();
            testCase.verifyWarning(@() gsua_dmatrix(T, 50, 'Method','Joint', ...
                'Pool', pool(:,1:3)), 'gsua_dmatrix:PoolTooSmall');
        end

        function tinyPoolFallsBackToMarginalSampling(testCase)
            [T,pool] = jointTestFixture();
            id = 'gsua_dmatrix:PoolTooSmall';
            warning('off', id);
            testCase.addTeardown(@() warning('on', id));

            M = gsua_dmatrix(T, 400, 'Method','Joint', 'Pool', pool(:,1:3), 'Seed', 4);
            testCase.verifySize(M, [400 size(T,1)]);
            % A marginal fallback spans the table range and does NOT reproduce the pool's
            % correlation -- the point is that it degrades honestly rather than pretending.
            testCase.verifyLessThan(abs(corr(M(:,1),M(:,2))), 0.2);
            testCase.verifyGreaterThanOrEqual(min(M(:,1)), T.Range(1,1));
            testCase.verifyLessThanOrEqual(max(M(:,1)), T.Range(1,2));
        end

        function adequatePoolDoesNotWarn(testCase)
            [T,pool] = jointTestFixture();
            testCase.verifyWarningFree(@() gsua_dmatrix(T, 50, 'Method','Joint', ...
                'Pool', pool(:,1:5), 'Seed', 4));
        end

        function minPoolNIsConfigurable(testCase)
            [T,pool] = jointTestFixture();
            testCase.verifyWarningFree(@() gsua_dmatrix(T, 50, 'Method','Joint', ...
                'Pool', pool(:,1:3), 'MinPoolN', 3, 'Seed', 4));
            testCase.verifyWarning(@() gsua_dmatrix(T, 50, 'Method','Joint', ...
                'Pool', pool(:,1:8), 'MinPoolN', 20), 'gsua_dmatrix:PoolTooSmall');
        end

        function nonFinitePoolColumnsAreDropped(testCase)
            [T,pool] = jointTestFixture();
            pool(2,4) = NaN;
            testCase.verifyWarning(@() gsua_dmatrix(T, 50, 'Method','Joint', ...
                'Pool', pool, 'Seed', 4), 'gsua_dmatrix:PoolNonFinite');

            id = 'gsua_dmatrix:PoolNonFinite';
            warning('off', id);
            testCase.addTeardown(@() warning('on', id));
            M = gsua_dmatrix(T, 200, 'Method','Joint', 'Pool', pool, 'Seed', 4);
            testCase.verifyTrue(all(isfinite(M(:))));
        end

        function allNonFinitePoolErrors(testCase)
            [T,pool] = jointTestFixture();
            pool(1,:) = NaN;
            testCase.verifyError(@() gsua_dmatrix(T, 50, 'Method','Joint', 'Pool', pool), ...
                'gsua_dmatrix:PoolAllNonFinite');
        end
    end
end
