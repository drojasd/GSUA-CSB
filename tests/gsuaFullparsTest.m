classdef gsuaFullparsTest < matlab.unittest.TestCase
    %GSUAFULLPARSTEST gsua_fullpars reinserts fixed factors into a free-parameter vector.

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
        function reinsertsFixedFactor(testCase)
            % Factor 3 (c) pinned at 0.05 -> dropped from T; gsua_fullpars must put it back in place.
            ranges = [500 2000; 0.01 0.5; 0.05 0.05; 0.01 1];
            [T,~] = gsua_dataprep('twoScaleUserModel', ranges, 'names',{'a','b','c','d'}, ...
                'out_names',{'count','frac'});
            testCase.verifyEqual(height(T),3, 'The pinned factor must be dropped from T.');

            [full,names] = gsua_fullpars(T,[1000; 0.1; 0.3]);   % free = a,b,d
            testCase.verifyEqual(full, [1000; 0.1; 0.05; 0.3], 'AbsTol',1e-12);
            testCase.verifyEqual(reshape(names,1,[]), {'a','b','c','d'});
        end

        function expandsAnEnsembleColumnwise(testCase)
            ranges = [500 2000; 0.01 0.5; 0.05 0.05; 0.01 1];
            [T,~] = gsua_dataprep('twoScaleUserModel', ranges, 'names',{'a','b','c','d'});
            free = [1 2 3; 10 20 30; 100 200 300];   % 3 free x 3 columns
            full = gsua_fullpars(T, free);
            testCase.verifyEqual(size(full),[4 3]);
            testCase.verifyEqual(full(3,:), [0.05 0.05 0.05], 'AbsTol',1e-12);  % fixed row
            testCase.verifyEqual(full([1 2 4],:), free);                         % free rows preserved
        end

        function wrongLengthErrors(testCase)
            ranges = [500 2000; 0.01 0.5; 0.05 0.05; 0.01 1];
            [T,~] = gsua_dataprep('twoScaleUserModel', ranges);
            testCase.verifyError(@() gsua_fullpars(T,[1;2]), 'GSUA:gsua_fullpars:BadLength');
        end
    end
end
