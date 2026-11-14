classdef gsuaShapleyTest < matlab.unittest.TestCase
    %GSUASHAPLEYTEST Shapley effects in gsua_sa ('SensMethod','Shapley'), both estimators.
    %
    %   Checked against analytic values:
    %   - Ishigami (independent uniform factors): Sh = [0.4358 0.4424 0.1218], between the
    %     first-order S = [0.3139 0.4424 0] and total ST = [0.5576 0.4424 0.2437] indices.
    %   - y = x1 + 2*x2 with corr(x1,x2) = 0.5 ('Given' with a correlated 'Pool'): Sh = [0.3393 0.6607].
    %   - Ishigami with x3 fixed at 1: Sh = [0.0899 0.9101].
    %   Both estimators telescope, so the indices sum to exactly 1.

    properties (Constant)
        ShIshigami = [0.4358 0.4424 0.1218];
        SIshigami  = [0.3139 0.4424 0];
        STIshigami = [0.5576 0.4424 0.2437];
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot,'Functions'), IncludingSubfolders=false));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                testDir, IncludingSubfolders=false));
        end
    end

    methods (Access = private)
        function T = ishigamiTable(~, ranges)
            if nargin < 2
                ranges = repmat([-pi pi],3,1);
            end
            [T,~] = gsua_dataprep('ishigamiUserModel', ranges, 'names',{'x1','x2','x3'});
        end
        function [Ts, Y, J] = runShapley(~, M, T, varargin)
            [Ts, Y, J] = gsua_sa(M, T, 'SensMethod','Shapley', 'parallel',false, 'verbose',false, ...
                varargin{:});
        end
    end

    methods (Test)
        function independentMatchesIshigami(testCase)
            T = testCase.ishigamiTable();
            rng(11); M = gsua_dmatrix(T, 8000);
            Ts = testCase.runShapley(M, T);
            testCase.verifyEqual(Ts.Sh_vec(:).', testCase.ShIshigami, 'AbsTol', 0.04);
            testCase.verifyEqual(sum(Ts.Sh_vec(:)), 1, 'AbsTol', 1e-10);
            testCase.verifyEqual(sum(Ts.Sh), 1, 'AbsTol', 1e-10);
        end

        function givenMatchesIshigami(testCase)
            T = testCase.ishigamiTable();
            rng(12); M = gsua_dmatrix(T, 4000);
            Ts = testCase.runShapley(M, T, 'ShapleyType','Given');
            testCase.verifyEqual(Ts.Sh_vec(:).', testCase.ShIshigami, 'AbsTol', 0.03);
            testCase.verifyEqual(sum(Ts.Sh_vec(:)), 1, 'AbsTol', 1e-10);
            testCase.verifyEqual(sum(Ts.Sh), 1, 'AbsTol', 1e-10);
        end

        function shapleyLiesBetweenFirstAndTotalOrder(testCase)
            % For independent factors Si <= Sh_i <= STi (up to sampling noise).
            T = testCase.ishigamiTable();
            rng(13); M = gsua_dmatrix(T, 4000);
            Sh = testCase.runShapley(M, T, 'ShapleyType','Given').Sh_vec(:).';
            testCase.verifyGreaterThanOrEqual(Sh + 0.03, testCase.SIshigami);
            testCase.verifyLessThanOrEqual(Sh - 0.03, testCase.STIshigami);
        end

        function givenUsesTheDependenceOfACorrelatedPool(testCase)
            [T,~] = gsua_dataprep('linearUserModel', [-4 4; -4 4], 'names',{'x1','x2'});
            rng(14); L = chol([1 .5; .5 1],'lower'); pool = L*randn(2,4000);   % Np x nPool
            Ts = testCase.runShapley(gsua_dmatrix(T,10), T, 'ShapleyType','Given', 'Pool',pool);
            testCase.verifyEqual(Ts.Sh_vec(:).', [0.3393 0.6607], 'AbsTol', 0.03);
            % The same model with independent factors shares the variance differently (1/5, 4/5):
            % the dependence, not only the model, decides the attribution.
            rng(15); Ti = testCase.runShapley(gsua_dmatrix(T,4000), T);
            testCase.verifyEqual(Ti.Sh_vec(:).', [0.2 0.8], 'AbsTol', 0.04);
        end

        function fixedFactorIsExcluded(testCase)
            T = testCase.ishigamiTable([-pi pi; -pi pi; 1 1]);
            rng(16); M = gsua_dmatrix(T, 4000);
            for type = {'Independent','Given'}
                Ts = testCase.runShapley(M, T, 'ShapleyType',type{1});
                testCase.verifyEqual(numel(Ts.Sh), height(T));
                testCase.verifyEqual(Ts.Sh_vec(:).', [0.0899 0.9101], 'AbsTol', 0.05, type{1});
            end
        end

        function duplicatedPoolVectorsAreEvaluatedOnce(testCase)
            T = testCase.ishigamiTable();
            rng(17); base = gsua_dmatrix(T, 300);
            pool = [base; base(1:100,:)].';                     % 100 repeated vectors
            [~, Y] = testCase.runShapley(base, T, 'ShapleyType','Given', 'Pool',pool);
            testCase.verifyEqual(size(Y,1), 300);
        end

        function methodNameIsCaseInsensitive(testCase)
            T = testCase.ishigamiTable();
            rng(18); M = gsua_dmatrix(T, 400);
            Ts = gsua_sa(M, T, 'SensMethod','shapley', 'parallel',false, 'verbose',false);
            testCase.verifyEqual(Ts.Properties.CustomProperties.SensMethod, 'Shapley');
            testCase.verifyTrue(all(isfinite(Ts.Sh)));
        end

        function badPoolShapeErrors(testCase)
            T = testCase.ishigamiTable();
            M = gsua_dmatrix(T, 20);
            testCase.verifyError(@() testCase.runShapley(M, T, 'ShapleyType','Given', ...
                'Pool',rand(5,7)), 'gsua_sa:ShapleyPoolSize');
        end

        function multiObjectiveIsRejected(testCase)
            T = testCase.ishigamiTable();
            M = gsua_dmatrix(T, 20);
            testCase.verifyError(@() testCase.runShapley(M, T, 'multiObjective',true), ...
                'gsua_sa:ShapleyMultiObjective');
        end
    end
end
