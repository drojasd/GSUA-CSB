classdef gsuaDataprepFlexibleNamesTest < matlab.unittest.TestCase
    %GSUADATAPREPFLEXIBLENAMESTEST names/out_names accept string arrays and char, not just cellstr.
    %
    %   gsua_userdefined (and sens_dataprep) validated 'names'/'out_names' with @iscellstr, which
    %   rejected the modern string-array form ["a" "b"] and a bare 'y' char. They now accept
    %   cellstr, string array, or char and normalize internally, so a caller is no longer forced
    %   into {'a','b'} form. Empty still auto-generates numeric names.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot, 'Functions'), IncludingSubfolders=false));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot, 'Examples'), IncludingSubfolders=false));
        end
    end

    properties
        R = [0.6 3.0; 0.05 0.5; 5 40];
    end

    methods (Test)
        function stringArrayNames(testCase)
            T = gsua_userdefined('pkAbsorptionModel', testCase.R, ...
                'names', ["ka","ke","V"], 'domain', [0 24]);
            testCase.verifyEqual(T.Properties.RowNames, {'ka';'ke';'V'});
        end

        function charOutName(testCase)
            T = gsua_userdefined('pkAbsorptionModel', testCase.R, ...
                'out_names', 'concentration', 'domain', [0 24]);
            testCase.verifyEqual(cellstr(T.Properties.CustomProperties.Vars), {'concentration'});
        end

        function cellstrStillWorks(testCase)
            % Back-compatibility: the previously-required form must still be accepted.
            T = gsua_userdefined('pkAbsorptionModel', testCase.R, ...
                'names', {'a','b','c'}, 'domain', [0 24]);
            testCase.verifyEqual(T.Properties.RowNames, {'a';'b';'c'});
        end

        function omittedNamesAutoGenerate(testCase)
            T = gsua_userdefined('pkAbsorptionModel', testCase.R, 'domain', [0 24]);
            testCase.verifyEqual(T.Properties.RowNames, {'1';'2';'3'});
        end

        function stringArrayNamesFlowThroughGsuaDataprep(testCase)
            % The dispatcher must forward the flexible form to gsua_userdefined unchanged.
            T = gsua_dataprep('pkAbsorptionModel', testCase.R, 'names', ["p1","p2","p3"], ...
                'domain', [0 24]);
            testCase.verifyEqual(T.Properties.RowNames, {'p1';'p2';'p3'});
        end
    end
end
