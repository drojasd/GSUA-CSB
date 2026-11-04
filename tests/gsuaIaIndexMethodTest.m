classdef gsuaIaIndexMethodTest < matlab.unittest.TestCase
    %GSUAIAINDEXMETHODTEST Identifiability index 'spread' vs legacy 'ci' (#4).
    %
    %   A parameter scattered across its whole prior but uncorrelated and with a well-pinned median
    %   must be flagged as poorly identified by the default 'spread' index, where the legacy 'ci'
    %   index (median-CI width, which shrinks with sample size) rates it as comparatively identified.
    %   gsua_ia and gsua_dia must agree on the shared index computation.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            testDir = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(toolboxRoot,'Functions'), IncludingSubfolders=false));
        end
    end

    methods (Test)
        function flatParameterFlaggedBySpreadNotCi(testCase)
            rng(42);
            a = 5 + 0.01*randn(1,300);   % tightly identified
            b = 10*rand(1,300);          % flat across the whole prior, uncorrelated
            Est = [a; b];
            T = table([0 10; 0 10],'RowNames',{'a','b'},'VariableNames',{'Range'});

            Ts = gsua_dia(T,Est,false,false,'IndexMethod','spread');
            Tc = gsua_dia(T,Est,false,false,'IndexMethod','ci');

            testCase.verifyGreaterThan(Ts.index(2),Tc.index(2));   % spread exposes the flat param
            testCase.verifyGreaterThan(Ts.index(2),0.5);           % flagged poorly identified
            testCase.verifyLessThan(Ts.index(1),0.4);              % tight param stays low
        end

        function iaMatchesDiaIndex(testCase)
            rng(7);
            Est = [5 + 0.01*randn(1,200); 10*rand(1,200)];
            T = table([0 10; 0 10],'RowNames',{'a','b'},'VariableNames',{'Range'});

            Td = gsua_dia(T,Est,false,false,'IndexMethod','spread');

            figs0 = findall(0,'Type','figure');
            testCase.addTeardown(@() close(setdiff(findall(0,'Type','figure'),figs0)));
            set(groot,'defaultFigureVisible','off');
            testCase.addTeardown(@() set(groot,'defaultFigureVisible','on'));
            Ti = gsua_ia(T,Est,false,false,false,false,'IndexMethod','spread');

            testCase.verifyEqual(Ti.index,Td.index,'AbsTol',1e-10);
        end
    end
end
