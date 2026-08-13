classdef gsuaNoisefloorNb2ParamTest < matlab.unittest.TestCase
    %GSUANOISEFLOORNB2PARAMTEST Verifies gsua_noisefloor's documented NB2 reparameterization.
    %
    %   gsua_noisefloor draws NB-distributed noise with a target mean mu and target dispersion k
    %   (var = mu + mu^2/k) via MATLAB's nbinrnd(R,P) with R=k, P=k/(k+mu) -- this is the exact
    %   formula documented in gsua_noisefloor's help text and used internally. Tested directly
    %   against nbinrnd (a public, real MATLAB function), not by reaching into gsua_noisefloor's
    %   private bootstrap-generation helper.

    properties (TestParameter)
        muK = struct( ...
            'moderate', struct('mu', 50, 'k', 8), ...
            'small_mu', struct('mu', 2, 'k', 8), ...
            'large_k', struct('mu', 100, 'k', 200));
    end

    methods (TestMethodSetup)
        function resetRandomSeed(testCase)
            originalRng = rng;
            testCase.addTeardown(@() rng(originalRng));
            rng(0, 'twister');
        end
    end

    methods (Test)
        function nb2MeanAndVarianceMatchTarget(testCase, muK)
            mu = muK.mu;
            k = muK.k;
            n = 200000;
            R = k;
            P = k / (k + mu);
            sample = nbinrnd(R, P, n, 1);

            targetVar = mu + mu^2 / k;
            testCase.verifyEqual(mean(sample), mu, 'RelTol', 0.02);
            testCase.verifyEqual(var(sample), targetVar, 'RelTol', 0.05);
        end
    end
end
