classdef gsuaClusterSummaryTest < matlab.unittest.TestCase
    %GSUACLUSTERSUMMARYTEST Per-cluster summary table + identifiability verdict (#8).

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot,'Functions'), IncludingSubfolders=false));
        end
    end

    methods (Test)
        function twoClusterSummary(testCase)
            rng(5);
            basin1 = [2+0.05*randn(60,1), 2+0.05*randn(60,1)];
            basin2 = [8+0.05*randn(60,1), 8+0.05*randn(60,1)];
            Est = [basin1; basin2]';
            T = table([0 10;0 10],'RowNames',{'a','b'},'VariableNames',{'Range'});
            [~,ci] = gsua_dia(T,Est,false,true,5,0.5);

            S = gsua_clustersummary(T,ci);
            testCase.verifyEqual(height(S),4);                      % 2 clusters x 2 params
            testCase.verifyTrue(all(strcmp(S.Verdict,'identifiable')));
            testCase.verifyTrue(all(S.MedianCI_low<=S.Median & S.Median<=S.MedianCI_high));
            testCase.verifyTrue(all(S.P2_5<=S.Median & S.Median<=S.P97_5));
        end

        function flatParameterNotIdentifiable(testCase)
            rng(1);
            Est = [5+0.01*randn(1,200); 10*rand(1,200)];
            T = table([0 10;0 10],'RowNames',{'a','b'},'VariableNames',{'Range'});
            [~,ci] = gsua_dia(T,Est,false,false);

            S = gsua_clustersummary(T,ci);
            a = S(strcmp(S.Parameter,'a'),:);
            b = S(strcmp(S.Parameter,'b'),:);
            testCase.verifyEqual(a.Verdict{1},'identifiable');
            testCase.verifyEqual(b.Verdict{1},'not identifiable');
            testCase.verifyGreaterThan(b.PctWidthRel,0.5);
        end
    end
end
