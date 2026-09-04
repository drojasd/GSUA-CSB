classdef gsuaDmatrixJointOrientationTest < matlab.unittest.TestCase
    %GSUADMATRIXJOINTORIENTATIONTEST Joint sampling must honor gsua_dmatrix's output contract.
    %
    %   Every existing gsua_dmatrix branch returns an N x Np matrix whose columns are
    %   row-aligned with the table and whose fixed-parameter columns are constant. The joint
    %   branch reads a Np x nPool pool and therefore has to transpose; these tests guard that
    %   it transposes the right way, still pins fixed parameters, and resolves its pool from
    %   the expected places.

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
        function outputIsNbyNp(testCase)
            T = jointTestFixture();
            M = gsua_dmatrix(T, 250, 'Method','Joint', 'Seed', 7);
            testCase.verifySize(M, [250 size(T,1)]);
            testCase.verifyTrue(all(isfinite(M(:))));
        end

        function everyBootstrapDrawIsAPoolColumn(testCase)
            % Bootstrap draws must be pool vectors verbatim -- drawing them intact is what
            % preserves EVERY dependence in the pool, not merely its pairwise correlations.
            [T,pool] = jointTestFixture();
            M = gsua_dmatrix(T, 200, 'Method','Joint', 'Seed', 3);
            worst = 0;
            for i = 1:size(M,1)
                worst = max(worst, min(sqrt(sum((pool - M(i,:)').^2, 1))));
            end
            testCase.verifyLessThan(worst, 1e-12);
        end

        function fixedParametersStayFixed(testCase)
            [T,pool] = jointTestFixture();
            T.Range(3,:) = [0.5 0.5];
            M = gsua_dmatrix(T, 100, 'Method','Joint', 'Pool', pool, 'Seed', 11);
            testCase.verifyEqual(M(:,3), repmat(0.5,100,1), 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(std(M(:,1)), 0);
        end

        function explicitPoolOverridesTableEst(testCase)
            T = jointTestFixture();
            alt = repmat([1.25; 1.6; 0.5], 1, 12);
            M = gsua_dmatrix(T, 50, 'Method','Joint', 'Pool', alt, 'Seed', 5);
            testCase.verifyEqual(M, repmat([1.25 1.6 0.5], 50, 1), 'AbsTol', 1e-12);
        end

        function wrongSizedPoolErrors(testCase)
            T = jointTestFixture();
            testCase.verifyError(@() gsua_dmatrix(T, 10, 'Method','Joint', ...
                'Pool', ones(2,20)), 'gsua_dmatrix:PoolSizeMismatch');
        end

        function missingPoolErrors(testCase)
            T = gsua_userdefined('jointManifoldFixtureFunc', [0.3 3; 0.3 3; 0.1 1.2]);
            testCase.verifyError(@() gsua_dmatrix(T, 10, 'Method','Joint'), ...
                'gsua_dmatrix:NoPool');
        end

        function jointTypeWithoutJointErrors(testCase)
            T = jointTestFixture();
            testCase.verifyError(@() gsua_dmatrix(T, 10, 'JointType','Gaussian'), ...
                'gsua_dmatrix:JointTypeWithoutJoint');
        end

        function marginalMethodsAreUnaffected(testCase)
            % Regression guard: adding a Method value must not disturb the existing ones.
            T = jointTestFixture();
            for m = {'LatinHypercube','Uniform','Sobol'}
                M = gsua_dmatrix(T, 40, 'Method', m{1});
                testCase.verifySize(M, [40 size(T,1)]);
                testCase.verifyTrue(all(M(:,1) >= T.Range(1,1) & M(:,1) <= T.Range(1,2)));
            end
        end
    end
end
