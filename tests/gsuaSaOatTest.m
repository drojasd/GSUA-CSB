classdef gsuaSaOatTest < matlab.unittest.TestCase
    %GSUASAOATTEST gsua_sa 'SensMethod','OAT' (one factor at a time).
    %
    %   Before v1.19 the single-objective OAT passed the single column M(:,i) to the model: every
    %   evaluation failed (left as zeros) and the call then crashed on an undefined STi_vec. OAT now
    %   varies factor i over its sampled values with the others at T.Nominal (the Python definition).
    %   For y = x1 + 2*x2 on [-4,4]^2 with nominal 0 the OAT cost of factor k is c_k^2 x_k^2, so
    %   Si = c^4/sum(c^4) = [1/17 16/17].

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

    methods (Test)
        function linearModelMatchesAnalyticShares(testCase)
            [T,~] = gsua_dataprep('linearUserModel', [-4 4; -4 4], 'names',{'x1','x2'});
            rng(21); M = gsua_dmatrix(T, 4000);
            [To, Y, J] = testCase.verifyWarningFree(@() gsua_sa(M, T, 'SensMethod','OAT', ...
                'parallel',false, 'verbose',false));
            testCase.verifyEqual(To.Si(:).', [1/17 16/17], 'AbsTol', 0.01);
            testCase.verifyEqual(sum(To.Si), 1, 'AbsTol', 1e-12);
            % STi is the OAT variance itself: check factor 1 by hand.
            p = repmat(T.Nominal(:).', size(M,1), 1); p(:,1) = M(:,1);
            cost = sum((linearUserModel(p) - linearUserModel(T.Nominal(:).')).^2, 2);
            testCase.verifyEqual(To.STi(1), var(cost,1), 'RelTol', 1e-10);
            % J is the cost of the base sample, one value per run (as for every method).
            testCase.verifySize(J, [size(Y,1) 1]);
            testCase.verifyFalse(any(strcmp(To.Properties.VariableNames,'Si_vec')), ...
                'OAT has no time-dependent indices.');
        end

        function timeDependentModelRunsWithoutFailedEvaluations(testCase)
            [T,~,~] = sirSencilloFixture();
            rng(22); M = gsua_dmatrix(T, 60);
            To = testCase.verifyWarningFree(@() gsua_sa(M, T, 'SensMethod','OAT', ...
                'parallel',false, 'verbose',false));
            testCase.verifyEqual(numel(To.Si), height(T));
            testCase.verifyTrue(all(isfinite(To.Si)) && all(To.Si >= 0));
            testCase.verifyEqual(sum(To.Si), 1, 'AbsTol', 1e-12);
        end
    end
end
