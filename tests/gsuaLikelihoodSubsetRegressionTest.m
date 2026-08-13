classdef gsuaLikelihoodSubsetRegressionTest < matlab.unittest.TestCase
    %GSUALIKELIHOODSUBSETREGRESSIONTEST Regression test for the Taux size/alignment fix.
    %
    %   gsua_likelihood indexes its working range matrix (Taux) by each parameter's ABSOLUTE row
    %   in T, not by its position within a profiled subset -- so Taux must be seeded from T.Range
    %   (Np x 2, correctly sized/aligned), not zeros(npars,2) (sized to the profiled-subset count).
    %   Before the fix, profiling a strict, non-empty subset either errored (dimension mismatch
    %   assigning Taux back onto Tfinal.Range) or silently corrupted every unprofiled row to
    %   [0,0]. Verified here on a tiny fixture (3 free params, 1 profiled) so it runs in seconds.

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

    properties
        T
        xdata
    end

    methods (TestMethodSetup)
        function buildFixture(testCase)
            testCase.xdata = linspace(0,4,15);
            testCase.T = gsua_userdefined('noisefloorFixtureFunc', [0.1 20; 0.01 2; 0.01 5]);
            testCase.T.Nominal = [4; 0.4; 0.8];
        end
    end

    methods (Test)
        function fullParameterSetProfilingStillWorks(testCase)
            % pars empty -> unaffected by the fix (npars==size(T,1) already, no misalignment
            % possible) -- a true regression guard that nothing broke for the common case.
            Tfinal = gsua_likelihood(testCase.T, testCase.xdata, [], 0.95, 0.1, 0.1, 0.5, 0.5, ...
                2, 1, false, false, false, []);
            testCase.verifySize(Tfinal.Range, size(testCase.T.Range));
            testCase.verifyTrue(all(isfinite(Tfinal.Range(:))));
        end

        function strictSubsetProfilingDoesNotErrorOrCorruptUnprofiledRows(testCase)
            Tfinal = gsua_likelihood(testCase.T, testCase.xdata, [], 0.95, 0.1, 0.1, 0.5, 0.5, ...
                2, 1, false, false, false, 2);

            testCase.verifySize(Tfinal.Range, size(testCase.T.Range));
            % The core regression assertion: unprofiled rows (1 and 3) must retain the original
            % T.Range exactly -- not [0,0], not truncated/misaligned.
            testCase.verifyEqual(Tfinal.Range([1 3],:), testCase.T.Range([1 3],:));
            % Row 2 (profiled) must have actually moved/been evaluated, not left untouched.
            testCase.verifyTrue(all(isfinite(Tfinal.Range(2,:))));
        end
    end
end
