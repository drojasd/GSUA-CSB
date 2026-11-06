classdef gsuaBootstrapPerRowNoiseTest < matlab.unittest.TestCase
    %GSUABOOTSTRAPPERROWNOISETEST Per-row noise in gsua_bootstrap + the zero-cost-scale guard (#9).
    %
    %   Applying a count noise model to a 0-1 fraction row can draw it to all zeros, which zeroes that
    %   row's cost scale and made fmincon fail "undefined at initial point" until the 100-error limit,
    %   never saying which row. The fix is twofold: gsua_bootstrap accepts a PER-ROW noise cellstr (so
    %   a fraction row can use 'none'/'gaussian'), and gsua_pe/gsua_costf guard a zero/undefined
    %   regulator by excluding that row with a named warning instead of failing.

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
        function [T,xdata,ydata,opt] = fixture()
            ranges = [500 2000; 0.01 0.5; 0.001 0.1; 0.01 1];
            [T,~] = gsua_dataprep('twoScaleUserModel', ranges, 'names',{'a','b','c','d'}, ...
                'out_names',{'count','frac'});
            xdata = 1:40;
            ydata = gsua_deval([1000;0.1;0.03;0.3],T,xdata);
            opt = optimoptions('fmincon','Display','off');
        end
    end

    methods (Test)
        function perRowNoiseCellstrRunsAndIsRecorded(testCase)
            [T,xdata,ydata,opt] = gsuaBootstrapPerRowNoiseTest.fixture();
            rng(2); [T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.3,'N',4, ...
                'opt',opt,'save',false,'timer',false);
            [E,~,out] = gsua_bootstrap(T,xdata,ydata,T.Estfmincon,res,'B',4,'N',3, ...
                'noise',{'quasipoisson','gaussian'},'opt',opt,'seed',1,'timer',false);
            testCase.verifySize(E,[4 4]);
            testCase.verifyClass(out.noise,'cell');
            testCase.verifyEqual(out.noise,{'quasipoisson','gaussian'});
        end

        function wrongLengthNoiseErrors(testCase)
            [T,xdata,ydata,opt] = gsuaBootstrapPerRowNoiseTest.fixture();
            Est0 = repmat(T.Nominal,1,3); res = [1 2 3];
            testCase.verifyError(@() gsua_bootstrap(T,xdata,ydata,Est0,res, ...
                'noise',{'poisson','none','none'},'B',3,'N',2,'margin',0.3,'timer',false), ...
                'gsua_bootstrap:BadNoise');
        end

        function zeroRegulatorRowIsGuardedNotFatal(testCase)
            % An all-zero data row would otherwise make fmincon undefined; gsua_pe must warn
            % ZeroRegulator (naming the row) and still converge the other row.
            [T,xdata,ydata,opt] = gsuaBootstrapPerRowNoiseTest.fixture();
            ydata(2,:) = 0;   % degenerate fraction row
            rng(2);
            f = @() gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.3,'N',2, ...
                'opt',opt,'save',false,'timer',false);
            testCase.verifyWarning(f,'GSUA:gsua_pe:ZeroRegulator');
            rng(2); [~,res] = f();
            testCase.verifyTrue(any(isfinite(res)), 'The healthy row must still be fittable.');
        end
    end
end
