classdef gsuaNoisefloorShowTest < matlab.unittest.TestCase
    %GSUANOISEFLOORSHOWTEST gsua_noisefloor 'Show' threshold-selection plots run headless (#9).

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
        function showRunsAndReturnsResult(testCase)
            [T, xdata, ydata] = noisefloorTestFixture();
            Est = T.Estfmincon;
            rng(0);
            res = sort(rand(size(Est,2),1)) * 1e-3;   % some below/above the floor, for the plot

            figs0 = findall(0,'Type','figure');
            testCase.addTeardown(@() close(setdiff(findall(0,'Type','figure'),figs0)));
            set(groot,'defaultFigureVisible','off');
            testCase.addTeardown(@() set(groot,'defaultFigureVisible','on'));
            warning('off','gsua_noisefloor:QuasiPoissonDegenerate');
            warning('off','gsua_noisefloor:UnderdispersedNB');
            warning('off','gsua_noisefloor:BestFitExceedsNoiseFloor');
            testCase.addTeardown(@() warning('on','all'));

            out = gsua_noisefloor(T, xdata, ydata, Est, res, 'B', 50, 'seed', 4, 'Show', true);
            testCase.verifyTrue(isfield(out,'accepted'));
            testCase.verifyEqual(numel(out.accepted), size(Est,2));
        end
    end
end
