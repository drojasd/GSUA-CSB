classdef gsuaMultiOutputCycleTest < matlab.unittest.TestCase
    %GSUAMULTIOUTPUTCYCLETEST Acceptance test: a matrix-output user model runs the full GSUA cycle.
    %
    %   A one-argument user function that returns a [nSeries x nSteps] matrix (Kind 6) must run
    %   gsua_ua -> gsua_sa -> gsua_pe(fmincon, margin 0.1) -> gsua_dia WITHOUT any scalar adapter,
    %   even with two series at very different scales and NaN in one of them. Before v1.18 gsua_pardeval
    %   crashed on such a model. These tests pin the shapes the cycle produces (3-D Y and per-series
    %   Si_vec) and that every stage completes; a single-output model is covered too, to confirm the
    %   2-D path is unchanged.

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
        function multiOutputFullCycle(testCase)
            ranges = [500 2000; 0.01 0.5; 0.001 0.1; 0.01 1];
            [T,~] = gsua_dataprep('twoScaleUserModel', ranges, ...
                'names',{'a','b','c','d'}, 'out_names',{'count','frac'});
            cp = T.Properties.CustomProperties;
            testCase.verifyEqual(cp.OutSize,[2 40]);
            testCase.verifyEqual(numel(cp.output),2);

            xdata = 1:40;
            truth = [1000; 0.1; 0.03; 0.3];
            yf = gsua_deval(truth,T,xdata);                 % 2 x 40
            rng(1);
            ydata = yf;
            ydata(1,:) = poissrnd(max(yf(1,:),0));          % count-like noise on row 1
            ydata(2,:) = yf(2,:) + 0.002*randn(1,40);       % small gaussian noise on row 2
            ydata(2,[5 15 25]) = NaN;                        % NaN in ONE series

            % --- gsua_ua: 3-D ensemble ---
            rng(2); M = gsua_dmatrix(T,60);
            Yua = gsua_ua(M,T,'plot',false,'parallel',false,'verbose',false);
            testCase.verifyEqual(size(Yua),[60 40 2], 'gsua_ua must return N x nSteps x nSeries.');

            % --- gsua_sa: per-series time-dependent indices ---
            rng(2); Tsa = gsua_sa(M,T,'SensMethod','Xiao','parallel',false,'verbose',false);
            testCase.verifyEqual(size(Tsa.Si_vec),[4 40 2], 'gsua_sa must return Np x Nd x nSeries.');
            testCase.verifyEqual(size(Tsa.STi_vec),[4 40 2]);

            % --- gsua_pe: fmincon, margin 0.1, no scalar adapter, NaN-tolerant ---
            opt = optimoptions('fmincon','Display','off');
            rng(3); [Tpe,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.1, ...
                'N',8,'opt',opt,'save',false,'timer',false);
            testCase.verifyTrue(any(isfinite(res)), 'At least one estimate must converge.');
            testCase.verifyEqual(size(Tpe.Estfmincon,1),4);

            % --- gsua_dia on the estimate ensemble ---
            [~,ci] = gsua_dia(Tpe, Tpe.Estfmincon, false, false, false, false);
            testCase.verifyGreaterThanOrEqual(ci.NumClusters,1);
            testCase.verifyTrue(isfield(ci,'OrigRange'));
        end

        function singleOutputCycleUnchanged(testCase)
            % A single observed series (Kind 5) still flows through the same cycle with 2-D indices.
            [T,xdata,ydata] = sirSencilloFixture();
            rng(2); M = gsua_dmatrix(T,60);
            Yua = gsua_ua(M,T,'plot',false,'parallel',false,'verbose',false);
            testCase.verifyEqual(size(Yua,1),60);
            rng(2); Tsa = gsua_sa(M,T,'SensMethod','Xiao','parallel',false,'verbose',false);
            testCase.verifyEqual(size(Tsa.Si_vec,1),height(T));
            testCase.verifyEqual(ndims(Tsa.Si_vec),2, 'Single series -> 2-D Si_vec (Np x Nd).');
            opt = optimoptions('fmincon','Display','off');
            rng(3); [Tpe,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.1, ...
                'N',6,'opt',opt,'save',false,'timer',false);
            testCase.verifyTrue(any(isfinite(res)));
            [~,ci] = gsua_dia(Tpe, Tpe.Estfmincon, false, false, false, false);
            testCase.verifyGreaterThanOrEqual(ci.NumClusters,1);
        end
    end
end
