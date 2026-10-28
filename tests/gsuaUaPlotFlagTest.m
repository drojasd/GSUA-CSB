classdef gsuaUaPlotFlagTest < matlab.unittest.TestCase
    %GSUAUAPLOTFLAGTEST gsua_ua(...,'plot',false) opens no figures and returns the same Y.
    %
    %   gsua_ua always opened figure(1) plus one Monte-Carlo-filtering figure per output,
    %   so batch runs calling it many times piled up dozens of windows (76+ reported from the
    %   EAFIT seminar). The 'plot' pair turns the figures off; Y must not depend on it.

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

    methods (Test)
        function plotFalseOpensNoFiguresAndKeepsY(testCase)
            T = gsua_dataprep('pkAbsorptionModel', [0.6 3.0; 0.05 0.5; 5 40], ...
                'domain', [0 24], 'out_names', {'c'});
            M = gsua_dmatrix(T, 20);
            xdata = [0.5 1 2 4 8 12 24];

            before = findall(groot, 'Type', 'figure');
            yQuiet = gsua_ua(M, T, 'xdata', xdata, 'parallel', false, ...
                'verbose', false, 'plot', false);
            after = findall(groot, 'Type', 'figure');
            testCase.verifyEqual(numel(after), numel(before), ...
                'gsua_ua(...,''plot'',false) must not open any figure.');

            % Same Y as the default (plotting) call; close what that call opens.
            yPlot = gsua_ua(M, T, 'xdata', xdata, 'parallel', false, 'verbose', false);
            opened = setdiff(findall(groot, 'Type', 'figure'), after);
            testCase.addTeardown(@() close(opened));
            testCase.verifyNotEmpty(opened, 'The default must still plot.');
            testCase.verifyEqual(yQuiet, yPlot);
        end
    end
end
