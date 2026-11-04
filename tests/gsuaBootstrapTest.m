classdef gsuaBootstrapTest < matlab.unittest.TestCase
    %GSUABOOTSTRAPTEST Parametric-bootstrap ensemble (gsua_bootstrap).
    %
    %   Refitting the same data from many starts collapses to one optimum (multistart best S0 lands
    %   near 2982.7 on the discrete-SIR fixture, truth 5000); the parametric bootstrap must instead
    %   reveal the wide S0 uncertainty the single best fit hides, and must be reproducible by seed.

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
        function ensembleShapeAndReproducibility(testCase)
            [T,xdata,ydata] = sirSencilloFixture();
            opt = optimoptions('fmincon','Display','off');
            rng(2);
            [T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.3,'N',10, ...
                'opt',opt,'save',false,'timer',false);
            E1 = gsua_bootstrap(T,xdata,ydata,T.Estfmincon,res,'B',6,'solver','fmincon', ...
                'N',4,'opt',opt,'seed',7,'timer',false);
            E2 = gsua_bootstrap(T,xdata,ydata,T.Estfmincon,res,'B',6,'solver','fmincon', ...
                'N',4,'opt',opt,'seed',7,'timer',false);
            testCase.verifySize(E1,[6 6]);
            testCase.verifyEqual(E1,E2,'A fixed seed must reproduce the ensemble exactly.');
        end

        function bootstrapRevealsS0NonIdentifiability(testCase)
            [T,xdata,ydata] = sirSencilloFixture();
            opt = optimoptions('fmincon','Display','off');
            rng(2);
            [T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.3,'N',10, ...
                'opt',opt,'save',false,'timer',false);
            % Multistart alone collapses to a single best S0 well below the truth of 5000.
            testCase.verifyLessThan(T.Estfmincon(1,1), 4000);

            [E,costs] = gsua_bootstrap(T,xdata,ydata,T.Estfmincon,res,'B',10,'solver','fmincon', ...
                'N',6,'opt',opt,'seed',1,'timer',false);
            testCase.verifyNumElements(costs,10);
            testCase.verifyTrue(all(isfinite(costs)));
            % The bootstrap must reveal the wide S0 uncertainty the single best fit hides.
            s0spread = max(E(1,:)) - min(E(1,:));
            testCase.verifyGreaterThan(s0spread,1500);
        end

        function badCumulativeErrors(testCase)
            [T,xdata,ydata] = sirSencilloFixture();
            Est0 = repmat(T.Nominal,1,3);
            res = [1 2 3];
            % ydata has one output row, so a 2-element cumulative flag is invalid.
            testCase.verifyError(@() gsua_bootstrap(T,xdata,ydata,Est0,res, ...
                'margin',0.3,'alpha',2,'cumulative',[false false],'B',3,'N',2,'timer',false), ...
                'gsua_bootstrap:BadCumulative');
        end
    end
end
