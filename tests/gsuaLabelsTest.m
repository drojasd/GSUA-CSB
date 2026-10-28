classdef gsuaLabelsTest < matlab.unittest.TestCase
    %GSUALABELSTEST Parameter names written in LaTeX render in the toolbox's figures.
    %
    %   Figures used the 'tex' interpreter only, so a name copied from a paper as
    %   '$\hat{\beta}_h$' showed its dollar signs and raw commands, with a "valid
    %   interpreter syntax" warning. gsua_labels switches a name set to 'latex' when any
    %   name looks like LaTeX, and leaves TeX names exactly as before.

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

    methods (TestMethodTeardown)
        function closeFigures(~)
            close all
        end
    end

    methods (Access = private)
        function T = latexTable(~, names)
            T = gsua_userdefined('pkAbsorptionModel', [0.6 3.0; 0.05 0.5; 5 40], ...
                'names', names, 'domain', [0 24], 'out_names', {'$c(t)$'});
        end

        function verifyNoInterpreterWarning(testCase, fcn)
            lastwarn('');
            fcn();
            drawnow
            [msg, ~] = lastwarn;
            testCase.verifyFalse(contains(msg, 'interpreter syntax'), msg);
        end
    end

    methods (Test)
        % ---- gsua_labels itself ----
        function texNamesAreUnchanged(testCase)
            in = {'\phi'; 'C'; '\mu_{w,l}'; 'beta_h'};
            [out, interp] = gsua_labels(in);
            testCase.verifyEqual(out, in);
            testCase.verifyEqual(interp, 'tex');
        end

        function dollarNamesSwitchTheSetToLatex(testCase)
            [out, interp] = gsua_labels({'$\phi$', 'C', '\mu_{w,l}'});
            testCase.verifyEqual(interp, 'latex');
            % already-delimited kept, TeX math wrapped, plain text left alone
            testCase.verifyEqual(out, {'$\phi$', 'C', '$\mu_{w,l}$'});
        end

        function latexOnlyCommandWithoutDollarsIsWrapped(testCase)
            [out, interp] = gsua_labels({'\hat{\beta}_h', 'k'});
            testCase.verifyEqual(interp, 'latex');
            testCase.verifyEqual(out, {'$\hat{\beta}_h$', 'k'});
        end

        function texCommandPrefixIsNotMistakenForLatex(testCase)
            % \beta starts with "\be", \textit-like prefixes must match whole words only.
            [~, interp] = gsua_labels({'\beta', '\theta_{dot}'});
            testCase.verifyEqual(interp, 'tex');
        end

        function charInCharOut(testCase)
            [out, interp] = gsua_labels('$x_0$');
            testCase.verifyClass(out, 'char');
            testCase.verifyEqual(interp, 'latex');
        end

        % ---- figures ----
        function scatterRendersLatexNames(testCase)
            T = testCase.latexTable({'$k_a$', '$\hat{k}_e$', '$\mathrm{V}$'});
            testCase.verifyNoInterpreterWarning(@() gsua_dmatrix(T, 60, 'Show', 'on'));
            ax = findall(gcf, 'Type', 'axes');
            yl = arrayfun(@(a) a.YLabel, ax);
            yl = yl(arrayfun(@(l) ~isempty(l.String), yl));
            testCase.verifyEqual(unique({yl.Interpreter}), {'latex'});
        end

        function longRowLabelStaysInsideTheFigure(testCase)
            T = testCase.latexTable({'$k_a$', '$\mathrm{elimination\ rate\ constant}_{e}$', '$V$'});
            gsua_dmatrix(T, 60, 'Show', 'on');
            drawnow
            ax = findall(gcf, 'Type', 'axes');
            left = inf;
            for a = ax'
                lb = a.YLabel;
                if isempty(lb.String) || lb.Rotation ~= 0, continue, end
                u = lb.Units; lb.Units = 'normalized'; e = lb.Extent; lb.Units = u;
                left = min(left, a.Position(1) + e(1)*a.Position(3));
            end
            testCase.verifyGreaterThanOrEqual(left, 0, ...
                'A row label runs past the left edge of the figure.');
        end

        function sensitivityPlotsRenderLatexNames(testCase)
            T = testCase.latexTable({'$k_a$', '$\hat{k}_e$', '$\mathrm{V}$'});
            S = [0.7; 0.25; 0.05];
            testCase.verifyNoInterpreterWarning(@() gsua_plot('Pie', T, S));
            figure
            testCase.verifyNoInterpreterWarning(@() gsua_plot('Bar', T, S));
        end

        function uncertaintyAndFilteringRenderLatexNames(testCase)
            T = testCase.latexTable({'$k_a$', '$\hat{k}_e$', '$\mathrm{V}$'});
            M = gsua_dmatrix(T, 30);
            testCase.verifyNoInterpreterWarning(@() gsua_ua(M, T, ...
                'xdata', [1 2 4 8], 'parallel', false, 'verbose', false));
        end
    end
end
