classdef gsuaPlotScatterNoSdoTest < matlab.unittest.TestCase
    %GSUAPLOTSCATTERNOSDOTEST The parameter scatter matrix no longer needs Simulink Design Optimization.
    %
    %   gsua_plot('ScatterParameter',...) (used by gsua_dmatrix(...,'Show','on')) called
    %   sdo.scatterPlot, which required the Simulink Design Optimization toolbox and rendered TeX
    %   parameter labels (e.g. \tau) poorly. It now uses base MATLAB's plotmatrix and keeps the
    %   real parameter names as TeX-interpreted axis labels.

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
        function sourceHasNoSdoCall(testCase)
            % Read the project's copy explicitly -- which('gsua_plot') can resolve to an
            % installed add-on copy that still carries the old sdo.scatterPlot call.
            testDir = fileparts(mfilename('fullpath'));
            src = fileread(fullfile(fileparts(testDir), 'Functions', 'gsua_plot.m'));
            testCase.verifyFalse(contains(src, 'sdo.scatterPlot'), ...
                'gsua_plot must not call sdo.scatterPlot (Simulink Design Optimization).');
            testCase.verifyTrue(contains(src, 'plotmatrix'), ...
                'gsua_plot ScatterParameter should use base-MATLAB plotmatrix.');
        end

        function scatterRunsAndLabelsAreTeX(testCase)
            T = gsua_userdefined('user_dependent', [0.1 0.3; 10 14], ...
                'domain', [0 25], 'names', ["\tau","\lambda"]);
            M = gsua_dmatrix(T, 200);
            fig = figure('Visible','off');
            testCase.addTeardown(@() close(fig));
            % Must not error and must produce a scatter matrix of axes.
            gsua_plot('ScatterParameter', T, 'LatinHypercube', M);
            ax = findall(fig, 'Type', 'axes');
            testCase.verifyGreaterThanOrEqual(numel(ax), 4);   % 2x2 scatter matrix
            labels = strings(0,1);
            for a = ax'
                if ~isempty(a.XLabel.String), labels(end+1) = string(a.XLabel.String); end %#ok<AGROW>
                if ~isempty(a.YLabel.String), labels(end+1) = string(a.YLabel.String); end %#ok<AGROW>
            end
            testCase.verifyTrue(any(contains(labels, '\tau')));
            testCase.verifyTrue(any(contains(labels, '\lambda')));
        end
    end
end
