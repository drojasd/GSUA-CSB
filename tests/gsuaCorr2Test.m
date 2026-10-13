classdef gsuaCorr2Test < matlab.unittest.TestCase
    %GSUACORR2TEST gsua_corr2 reproduces corr2 without the Image Processing Toolbox.
    %
    %   gsua_corr2 is the toolbox-free 2-D Pearson correlation used by the cost and
    %   sensitivity functions (gsua_costf/costfMulti, gsua_covmetric, gsua_oatr(2),
    %   gsua_rcostf, gsua_corr2omitnan) in place of Image Processing Toolbox corr2.
    %   These checks pin its behavior to analytic values, so the test does not itself
    %   depend on corr2 being installed.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot, 'Functions'), IncludingSubfolders=false));
        end
    end

    methods (Test)
        function perfectPositiveIsOne(testCase)
            testCase.verifyEqual(gsua_corr2([1 2 3 4],[2 4 6 8]), 1, 'AbsTol', 1e-12);
        end

        function perfectNegativeIsMinusOne(testCase)
            testCase.verifyEqual(gsua_corr2([1 2 3 4],[8 6 4 2]), -1, 'AbsTol', 1e-12);
        end

        function constantInputIsNaN(testCase)
            % corr2 returns 0/0 = NaN when an input has zero variance; downstream
            % callers guard on isnan(...), so this must stay NaN, not 0.
            testCase.verifyTrue(isnan(gsua_corr2([5 5 5],[1 2 3])));
        end

        function operatesOverAllElementsLikeCorr2(testCase)
            % corr2 linearizes: a matrix must give the same answer as its vectorization.
            A = magic(4); B = reshape(1:16,4,4) + 0.3*reshape(sin(1:16),4,4);
            testCase.verifyEqual(gsua_corr2(A,B), gsua_corr2(A(:).',B(:).'), 'AbsTol', 1e-12);
        end

        function matchesHandComputedValue(testCase)
            a = [1 2 3 5];  b = [2 1 4 3];
            am = a - mean(a); bm = b - mean(b);
            expected = sum(am.*bm)/sqrt(sum(am.^2)*sum(bm.^2));
            testCase.verifyEqual(gsua_corr2(a,b), expected, 'AbsTol', 1e-12);
        end
    end
end
