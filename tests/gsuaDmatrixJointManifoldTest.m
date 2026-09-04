classdef gsuaDmatrixJointManifoldTest < matlab.unittest.TestCase
    %GSUADMATRIXJOINTMANIFOLDTEST The analytic control: does a sampler stay on the manifold?
    %
    %   jointTestFixture's pool lies exactly on the hyperbola a*b = 2, the manifold that
    %   y = a*b*exp(-k*t) actually identifies. Because that manifold is known in closed form,
    %   "did the sampler leave it" is a measurement, not a judgement call.
    %
    %   This is the test that justifies Bootstrap being the default JointType. The manifold is
    %   CURVED, and a Gaussian fitted to a curved ridge must place mass beside it -- so
    %   preserving pairwise linear correlation (which the Gaussian mode does well, see
    %   gsuaDmatrixJointCorrelationTest) is demonstrably NOT sufficient to stay on it. If
    %   gaussianLeavesTheManifold ever stops holding, the default should be revisited.

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

    methods (Static, Access = private)
        function leak = manifoldLeak(M)
            % Distance from the identified manifold a*b = 2.
            leak = abs(M(:,1).*M(:,2) - 2);
        end
    end

    methods (Test)
        function bootstrapStaysOnTheManifoldExactly(testCase)
            T = jointTestFixture();
            M = gsua_dmatrix(T, 1000, 'Method','Joint', 'Seed', 13);
            testCase.verifyLessThan( ...
                max(gsuaDmatrixJointManifoldTest.manifoldLeak(M)), 1e-10);
        end

        function gaussianLeavesTheManifold(testCase)
            T = jointTestFixture();
            Mg = gsua_dmatrix(T, 1000, 'Method','Joint', 'JointType','Gaussian', 'Seed', 13);
            Mb = gsua_dmatrix(T, 1000, 'Method','Joint', 'Seed', 13);
            leakG = mean(gsuaDmatrixJointManifoldTest.manifoldLeak(Mg));
            leakB = mean(gsuaDmatrixJointManifoldTest.manifoldLeak(Mb));
            testCase.verifyGreaterThan(leakG, 1e-3);
            testCase.verifyGreaterThan(leakG, 1e3*max(leakB, eps));
        end

        function marginalSamplingLeaksMostOfAll(testCase)
            T = jointTestFixture();
            Ml = gsua_dmatrix(T, 1000, 'Seed', 13);
            Mg = gsua_dmatrix(T, 1000, 'Method','Joint', 'JointType','Gaussian', 'Seed', 13);
            leakL = mean(gsuaDmatrixJointManifoldTest.manifoldLeak(Ml));
            leakG = mean(gsuaDmatrixJointManifoldTest.manifoldLeak(Mg));
            testCase.verifyGreaterThan(leakL, leakG);
        end

        function smoothBootstrapStaysCloserThanGaussian(testCase)
            % Smoothing perturbs draws off the manifold on purpose (that is how it produces
            % new points), but anchoring each draw to a real pool vector should keep it
            % closer than a single global Gaussian does.
            T = jointTestFixture();
            Ms = gsua_dmatrix(T, 1000, 'Method','Joint', ...
                'JointType','SmoothBootstrap', 'Seed', 13);
            Mg = gsua_dmatrix(T, 1000, 'Method','Joint', 'JointType','Gaussian', 'Seed', 13);
            testCase.verifyLessThan( ...
                mean(gsuaDmatrixJointManifoldTest.manifoldLeak(Ms)), ...
                mean(gsuaDmatrixJointManifoldTest.manifoldLeak(Mg)));
        end

        function smoothBootstrapReproducesPoolSpread(testCase)
            % The variance correction exists so smoothing reproduces the ensemble covariance
            % rather than inflating it by (1+h^2). Without it this ratio drifts above 1.
            [T,pool] = jointTestFixture();
            Ms = gsua_dmatrix(T, 20000, 'Method','Joint', ...
                'JointType','SmoothBootstrap', 'Seed', 99);
            testCase.verifyEqual(std(Ms(:,1))/std(pool(1,:)), 1, 'AbsTol', 0.06);
        end
    end
end
