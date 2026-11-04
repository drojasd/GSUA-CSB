classdef gsuaUaPlotOrientationTest < matlab.unittest.TestCase
    %GSUAUAPLOTORIENTATIONTEST gsua_plot('UncertaintyAnalysis') orientation (v1.17.1 fix).
    %
    %   plot(t,M) is orientation-ambiguous when M is square: MATLAB then plots M's COLUMNS. With the
    %   uncertainty ensemble laid out as N simulations x numel(t) time points, a square ensemble
    %   (N == numel(xdata)) therefore came out transposed -- each drawn "curve" was one time instant
    %   sweeping the simulations, not one simulation over time. These tests build the square case and
    %   verify every plotted line is a simulation. They also cover the column-xdata case (which must
    %   still draw curves, not a single-instant boxplot) and the nominal-curve handle.

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

    methods (Static)
        function [T,Y,t,len,N] = squareCase()
            % A valid GSUA table (for its CustomProperties) + a SQUARE ensemble whose row i is the
            % constant i, so a correct plot gives N constant lines and the transposed bug gives
            % increasing lines.
            T = sirSencilloFixture();
            len = 15; N = len;                 % square: N == numel(t)
            t = 1:len;
            Y = zeros(N,len,1);
            for i = 1:N
                Y(i,:,1) = i;                  % simulation i is the constant i over time
            end
        end
    end

    methods (Test)
        function squareEnsembleEachLineIsASimulation(testCase)
            [T,Y,t,len,N] = gsuaUaPlotOrientationTest.squareCase();
            fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            gsua_plot('UncertaintyAnalysis', T, Y, t);
            ax = findobj(fig,'Type','axes');
            lines = findobj(ax(1),'Type','line');
            testCase.verifyNumElements(lines, N, ...
                'A square ensemble must draw one line per simulation, not one per time instant.');
            vals = zeros(numel(lines),1);
            for i = 1:numel(lines)
                yd = lines(i).YData;
                testCase.verifyEqual(numel(yd), len, 'Each line must span the time axis.');
                testCase.verifyEqual(max(yd)-min(yd), 0, 'AbsTol', 1e-12, ...
                    'Each line must be one simulation (constant row), not a transposed time column.');
                vals(i) = yd(1);
            end
            testCase.verifyEqual(sort(vals), (1:N)', ...
                'Every simulation must appear exactly once.');
        end

        function columnXdataStillDrawsCurves(testCase)
            % A column xdata (numel > 1) must still select curves, not the single-instant boxplot.
            [T,Y,~,len,N] = gsuaUaPlotOrientationTest.squareCase();
            tcol = (1:len).';
            fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            gsua_plot('UncertaintyAnalysis', T, Y, tcol);
            ax = findobj(fig,'Type','axes');
            testCase.verifyEmpty(findobj(ax(1),'Tag','Box'), ...
                'A column xdata must not fall through to a boxplot.');
            lines = findobj(ax(1),'Type','line');
            testCase.verifyNumElements(lines, N);
            for i = 1:numel(lines)
                yd = lines(i).YData;
                testCase.verifyEqual(numel(yd), len);
                testCase.verifyEqual(max(yd)-min(yd), 0, 'AbsTol', 1e-12);
            end
        end

        function nominalHandleIsTheNominalCurve(testCase)
            % case 5: with a nominal overlay, handle N+1 must be the (red) nominal curve.
            [T,Y,t,len,N] = gsuaUaPlotOrientationTest.squareCase();
            yNom = 100*ones(1,len);            % distinct from every simulation row (1..N)
            fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            gsua_plot('UncertaintyAnalysis', T, Y, t, yNom);
            ax = findobj(fig,'Type','axes');
            lines = findobj(ax(1),'Type','line');
            testCase.verifyNumElements(lines, N+1, 'N simulations plus the nominal curve.');
            red = findobj(ax(1),'Type','line','Color',[1 0 0]);
            testCase.verifyNumElements(red, 1, 'Exactly one nominal (red) curve.');
            testCase.verifyEqual(red.YData, yNom, 'AbsTol', 1e-12, ...
                'The nominal handle must carry the nominal output, not a simulation.');
        end
    end
end
