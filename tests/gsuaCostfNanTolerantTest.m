classdef gsuaCostfNanTolerantTest < matlab.unittest.TestCase
    %GSUACOSTFNANTOLERANTTEST Regression tests for NaN-tolerant gsua_costf/gsua_pe.
    %
    %   Covers the feature request: gsua_pe's fmincon-family cost path
    %   (gsua_costf.m, and the regulator calculation at gsua_pe.m ~line 151)
    %   must ignore per-row NaN in ydata rather than let a single NaN
    %   poison the whole scalar cost -- needed for multi-output joint fits
    %   where different signals have independently gappy real-world data.

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
        Xdata
        Ydata       % row 1 has a NaN block, row 2 is fully populated
        Yfunction
        Margin = 1.1
        Alpha = 2
    end

    methods (TestMethodSetup)
        function buildFixtureData(testCase)
            testCase.Xdata = 1:10;
            testCase.Ydata = [ ...
                1.0, 2.0, NaN, NaN, NaN, 6.0, 7.0, 8.0, 9.0, 10.0; ...
                2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0];
            testCase.Yfunction = testCase.Ydata + 0.5;
            % Candidate model output at the positions where ydata is NaN is
            % excluded from every calculation regardless of its value (the
            % mask is driven by ydata's NaN pattern) -- give it a finite
            % stand-in anyway so a test failure can't be masked by NaN
            % arithmetic quietly producing NaN==NaN-is-false elsewhere.
            testCase.Yfunction(isnan(testCase.Ydata)) = 0;
        end
    end

    methods (Test)
        function noNanPropagatesToScalarCost(testCase)
            % (a) A NaN anywhere in ydata must not turn the final scalar
            % cost NaN -- it should only exclude that sample from its own
            % row's MSE/correlation terms.
            len = length(testCase.Xdata);
            inputs = size(testCase.Ydata, 1);
            regulator = gsuaCostfNanTolerantTest.referenceRegulator(testCase.Ydata, testCase.Margin);

            cost = gsua_costf(inputs, regulator, len, testCase.Ydata, testCase.Yfunction, testCase.Alpha);

            testCase.verifyFalse(isnan(cost));
            testCase.verifyTrue(isfinite(cost));
        end

        function nanRowCostMatchesManualNonNanRestriction(testCase)
            % (b) A row containing a NaN block must produce the same cost
            % as calling gsua_costf on that row alone, pre-restricted by
            % hand to its non-NaN entries -- i.e. the NaN-aware path and a
            % manual non-NaN restriction agree.
            row = 1;
            mask = ~isnan(testCase.Ydata(row, :));
            yRestricted = testCase.Ydata(row, mask);
            yfRestricted = testCase.Yfunction(row, mask);
            lenRestricted = numel(yRestricted);
            regulatorRestricted = gsuaCostfNanTolerantTest.referenceRegulator(yRestricted, testCase.Margin);

            costFromRestriction = gsua_costf(1, regulatorRestricted, lenRestricted, ...
                yRestricted, yfRestricted, testCase.Alpha);

            len = length(testCase.Xdata);
            regulatorFull = gsuaCostfNanTolerantTest.referenceRegulator(testCase.Ydata(row, :), testCase.Margin);
            costFromNanRow = gsua_costf(1, regulatorFull, len, ...
                testCase.Ydata(row, :), testCase.Yfunction(row, :), testCase.Alpha);

            testCase.verifyEqual(costFromNanRow, costFromRestriction, AbsTol=1e-10);
        end

        function nanFreeCallIsByteIdenticalToPreFixFormula(testCase)
            % (c) With no NaN present, results must be byte-identical to
            % the pre-fix formula (plain sum/len divisor, plain corr2) --
            % not just numerically close.
            yClean = testCase.Ydata;
            yClean(1, 3:5) = [3.0, 4.0, 5.0];
            yfClean = yClean + 0.5;
            len = length(testCase.Xdata);
            inputs = size(yClean, 1);

            regulatorNew = sum((yClean - yClean * testCase.Margin) .^ 2, 2, 'omitnan') ./ sum(~isnan(yClean), 2);
            regulatorOld = sum((yClean - yClean * testCase.Margin) .^ 2, 2) / len;
            testCase.verifyEqual(regulatorNew, regulatorOld);

            costNew = gsua_costf(inputs, regulatorNew, len, yClean, yfClean, testCase.Alpha);
            costOld = gsuaCostfNanTolerantTest.preFixCostf(inputs, regulatorOld, len, yClean, yfClean, testCase.Alpha);
            testCase.verifyEqual(costNew, costOld);
        end

        function corr2OmitNanMatchesCorr2WhenNoNan(testCase)
            % (d1) No NaN present -- must match corr2 exactly.
            a = [1.0, 2.0, 3.0, 4.0, 5.0];
            b = [2.0, 3.5, 5.0, 6.5, 8.0];
            testCase.verifyEqual(gsua_corr2omitnan(a, b), corr2(a, b));
        end

        function corr2OmitNanReturnsZeroBelowTwoValidPoints(testCase)
            % (d2) Fewer than 2 jointly-valid entries -- degrades to 0
            % (no correlation credit) instead of erroring or returning NaN.
            a = [1.0, NaN, NaN];
            b = [2.0, 3.0, 4.0];
            testCase.verifyEqual(gsua_corr2omitnan(a, b), 0);
        end

        function gsuaPeFminconCompletesWithNanInOneOutputRow(testCase)
            % (e) End-to-end: gsua_pe with solver='fmincon' and margin
            % activating the gsua_costf multi-objective path must not
            % crash with "Objective function is undefined at initial
            % point" when one output row has NaN gaps.
            truePars = [3, 0.5, 1.2];
            xdata = linspace(0, 5, 20); %#ok<NASGU> % domain is hardcoded inside nanPeE2eModelFunc to match
            ydata = nanPeE2eModelFunc(truePars);
            ydata(1, 6:10) = NaN;

            table = gsua_userdefined('nanPeE2eModelFunc', [0 10; 0 5; 0 5]);
            [~, res] = gsua_pe(table, linspace(0, 5, 20), ydata, 'N', 1, ...
                'solver', 'fmincon', 'margin', 0.1, 'ipoint', [2, 0.4, 1]);

            testCase.verifyFalse(isnan(res));
            testCase.verifyTrue(isfinite(res));
        end
    end

    methods (Static, Access = private)
        function regulator = referenceRegulator(ydata, margin)
            regulator = sum((ydata - ydata * margin) .^ 2, 2, 'omitnan') ./ sum(~isnan(ydata), 2);
        end

        function cos = preFixCostf(inputs, regulator, len, ydata, yfunction, alpha)
            % Reimplementation of gsua_costf's pre-fix formula (plain
            % sum/len, plain corr2), for byte-identical-on-clean-data
            % comparison against the current NaN-tolerant implementation.
            cost = (sum((ydata - yfunction) .^ 2, 2) / len) ./ regulator;
            for i = 1:inputs
                cost(i) = ((2 - corr2(ydata(i, :), yfunction(i, :))) * cost(i)) ^ alpha;
            end
            cos = sum(cost) / inputs;
        end
    end
end
