classdef gsuaIntrpOrientationTest < matlab.unittest.TestCase
    %GSUAINTRPORIENTATIONTEST gsua_intrp accepts any ode45 output orientation.
    %
    %   A user-defined model (Kind 5) can return its ODE solution in either of
    %   the two orientations ode45 produces -- sol = ode45(...) gives nStates x N,
    %   while [t,y] = ode45(...) gives N x nStates -- transposed or not, with or
    %   without sol.solver set. gsua_intrp now normalizes all of these to the
    %   canonical sol.x = 1xN / sol.y = nStates x N / sol.solver = '' form, so
    %   gsua_eval, gsua_pardeval and gsua_deval (the extra==1 and extra==2 paths)
    %   all return the same correct result regardless of how the model was written.

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
        nominal = [100; 10; 1; 0.5; 0.3; 0.1];   % [y0(1:3); k(1:3)]
        domain  = [0 10];
        xdata   = linspace(0, 10, 25);
    end

    properties (TestParameter)
        form = {'canonical','notsolver','ty','structTransposed'};
    end

    methods (Access = private)
        function T = buildTable(testCase, form)
            ranges = [0.9*testCase.nominal, 1.1*testCase.nominal];
            cfg = struct('form', form);
            T = gsua_dataprep('orientationFixtureModel', ranges, ...
                'domain', testCase.domain, 'nominal', testCase.nominal, ...
                'out_names', {'y1','y2','y3'}, 'opt', cfg);
        end

        function E = analytic(testCase)
            % Exact solution y_j(t) = y0_j * exp(-k_j * t), returned nStates x numel(xdata).
            y0 = testCase.nominal(1:3);
            k  = testCase.nominal(4:6);
            E  = y0 .* exp(-k .* testCase.xdata);
        end
    end

    methods (Test)
        function evalMatchesAnalytic(testCase, form)
            % gsua_eval path (gsua_intrp extra==1).
            T = testCase.buildTable(form);
            y = gsua_eval(testCase.nominal, T, testCase.xdata, [], false, false, false);
            % nSim==1, 3 outputs -> [1 x numel x 3]; drop the sim dimension.
            got = permute(squeeze(y), [2 1]);   % -> 3 x numel
            testCase.verifyEqual(got, testCase.analytic, 'RelTol', 5e-3, 'AbsTol', 1e-3);
        end

        function pardevalMatchesAnalytic(testCase, form)
            % gsua_pardeval batch path (gsua_intrp extra==1), 2 identical rows.
            T = testCase.buildTable(form);
            M = [testCase.nominal.'; testCase.nominal.'];
            y = gsua_pardeval(M, T, testCase.xdata, false);
            testCase.verifySize(y, [2, numel(testCase.xdata), 3]);
            got = permute(squeeze(y(1,:,:)), [2 1]);   % -> 3 x numel
            testCase.verifyEqual(got, testCase.analytic, 'RelTol', 5e-3, 'AbsTol', 1e-3);
        end

        function devalMatchesAnalytic(testCase, form)
            % gsua_deval single-point path (gsua_intrp extra==2), used by gsua_pe/likelihood.
            T = testCase.buildTable(form);
            y = gsua_deval(testCase.nominal.', T, testCase.xdata);
            testCase.verifySize(y, [3, numel(testCase.xdata)]);
            testCase.verifyEqual(y, testCase.analytic, 'RelTol', 5e-3, 'AbsTol', 1e-3);
        end
    end

    methods (Test)   % direct gsua_intrp unit checks
        function trueMismatchGivesClearError(testCase)
            % A sol.y whose neither dimension matches numel(sol.x) is a real error,
            % and must be reported clearly rather than as "X and V must be the same length".
            bad.x = linspace(0,1,20);
            bad.y = zeros(2, 25);      % 25 ~= 20 on both axes
            testCase.verifyError(@() gsua_intrp(bad, bad.x, 1:2), 'GSUA:gsua_intrp:BadShape');
        end

        function singleOutputVectorAccepted(testCase)
            % A one-state model returning a plain row/column vector is oriented to 1xN.
            s.x = linspace(0,1,15);
            s.y = (1:15).';           % 15x1 column
            out = gsua_intrp(s, s.x, 1);
            testCase.verifyEqual(out(:).', 1:15, 'RelTol', 1e-9);
        end
    end
end
