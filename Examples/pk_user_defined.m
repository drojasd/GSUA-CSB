%[text] # End-to-end system identification: a user-defined model
%[text] This example takes a pharmacokinetic model all the way from a plain MATLAB function to confidence intervals on its parameters, using a user-defined (Kind 6) model rather than Simulink or Symbolic Math.
%[text] The point of the example is not that the fit succeeds. It is that **a good fit tells you almost nothing about whether your parameters are identifiable**, and that the toolbox will tell you the difference if you ask it.
%[text] The model is a one-compartment pharmacokinetic model with first-order absorption, the standard description of a drug taken orally:
%[text] $ c(t)=\\frac{D\\,k_a}{V(k_a-k_e)}\\left(e^{-k_e t}-e^{-k_a t}\\right) $
%[text] with three factors to identify:
%[text] - $ k_a $ — absorption rate constant (1/h)
%[text] - $ k_e $ — elimination rate constant (1/h)
%[text] - $ V $ — apparent volume of distribution (L) \
%[text] The dose $ D=100 $ mg is known. The implementation lives in |Examples/pkAbsorptionModel.m|; it returns an ODE-solver-shaped struct so that |gsua_eval| can interpolate it onto any sampling times.
rng(0,'twister')
type pkAbsorptionModel.m
%%
%[text] ## 1. Preparing the environment
%[text] |gsua_dataprep| builds the summary table |T| that every other toolbox function consumes. For a user-defined model it needs the function name, the factor ranges, and the time domain.
ranges = [0.6 3.0;     % ka
          0.05 0.5;    % ke
          5    40];    % V
[T,~] = gsua_dataprep('pkAbsorptionModel', ranges, 'domain', [0 24], ...
    'names', {'ka','ke','V'}, 'out_names', {'concentration'});
T
%[text] The $ k_a $ range is deliberately kept above the $ k_e $ range. At $ k_a=k_e $ the closed form is singular, and swapping the two leaves $ c(t) $ unchanged — the classic *flip-flop* ambiguity. Excluding it keeps this example about experimental design rather than about an algebraic accident.
%%
%[text] ## 2. Synthetic data
%[text] Using synthetic data means the truth is known, so the confidence intervals can be checked rather than merely reported. The true factor values are $ k_a=1.2 $, $ k_e=0.25 $, $ V=15 $, and the measurements carry 8% proportional noise on a realistic sampling schedule.
truth = [1.2; 0.25; 15];
T.Nominal = truth;
xdata = [0.25 0.5 1 1.5 2 3 4 6 8 10 12 16 20 24];
clean = gsua_eval(T.Nominal, T, xdata);
ydata = clean .* (1 + 0.08*randn(size(clean)));
tdense = linspace(0.05, 24, 300);
plot(tdense, gsua_eval(truth, T, tdense), 'LineWidth', 1.5)
hold on
plot(xdata, ydata, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
hold off
xlabel('time (h)')
ylabel('concentration (mg/L)')
legend('true model','measurements','Location','northeast')
title('Simulated single-dose concentration data')
grid on
%%
%[text] ## 3. First attempt: estimate all three factors
%[text] |gsua_pe| runs a multistart estimation. |'margin',0.1| selects the correlation-penalized cost, which also records the margin on the table so downstream functions can recover it.
[T3,res3] = gsua_pe(T, xdata, ydata, 'solver','lsqc', 'N',20, 'margin',0.1, 'timer',false);
Est3 = T3.Estlsqc;
estimates = table(truth, Est3(:,1), 'VariableNames', {'true','estimated'}, ...
    'RowNames', T3.Properties.RowNames)
%[text] The fit is good and every one of the 20 multistart runs converged to the same cost, which is usually taken as a sign of a healthy estimation problem.
plot(tdense, gsua_eval(Est3(:,1), T3, tdense), 'LineWidth', 1.5)
hold on
plot(xdata, ydata, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
hold off
xlabel('time (h)')
ylabel('concentration (mg/L)')
legend('fitted model','measurements','Location','northeast')
title(sprintf('Fit with all three factors free (cost = %.4g)', min(res3)))
grid on
%%
%[text] ## 4. Diagnosis: is that fit identifiable?
%[text] The correlation between the repeated estimates is the first warning sign. Values near $ \\pm 1 $ mean the factors are trading off against each other: many different combinations reproduce the same curve.
corr3 = array2table(corr(Est3'), 'VariableNames', T3.Properties.RowNames, ...
    'RowNames', T3.Properties.RowNames)
%[text] |gsua_likelihood| turns that into a number you can act on. It profiles each factor — stepping it away from the estimate while re-fitting the others — and reports the interval where the fit remains statistically acceptable.
Tci3 = gsua_likelihood(T3, xdata, ydata, 0.95, 0.05, 0.1, 0.01, 0.01, 15, 1, false, false, false, []);
ci3 = table(Tci3.Range(:,1), Tci3.Range(:,2), Tci3.Range(:,2)-Tci3.Range(:,1), ranges(:,1), ranges(:,2), ...
    'VariableNames', {'CI_low','CI_high','width','prior_low','prior_high'}, ...
    'RowNames', T3.Properties.RowNames)
%[text] This is the result worth stopping on. The confidence interval for $ k_a $ spans its **entire prior range** — the data constrain it not at all — and $ k_e $ runs to its upper bound. A model that fits this well is still telling us that these three factors cannot be separated from a single oral concentration curve. That is a textbook result in pharmacokinetics, not a failure of the optimizer.
%%
%[text] ## 5. Remedy: fix what an independent experiment already knows
%[text] The usual resolution is to determine $ V $ separately, from an intravenous study where it is directly identifiable, and then estimate only the absorption and elimination constants. In the toolbox a factor is fixed by giving it a degenerate range; fixed factors drop out of |T| entirely.
rangesFixed = [0.6 3.0; 0.05 0.5; 15 15];
[Tf,~] = gsua_dataprep('pkAbsorptionModel', rangesFixed, 'domain', [0 24], ...
    'names', {'ka','ke','V'}, 'out_names', {'concentration'});
Tf.Nominal = truth(1:height(Tf));
[T2,res2] = gsua_pe(Tf, xdata, ydata, 'solver','lsqc', 'N',20, 'margin',0.1, 'timer',false);
Tci2 = gsua_likelihood(T2, xdata, ydata, 0.95, 0.05, 0.1, 0.01, 0.01, 15, 1, false, false, false, []);
ci2 = table(truth(1:2), T2.Estlsqc(:,1), Tci2.Range(:,1), Tci2.Range(:,2), ...
    Tci2.Range(:,2)-Tci2.Range(:,1), ...
    'VariableNames', {'true','estimated','CI_low','CI_high','width'}, ...
    'RowNames', T2.Properties.RowNames)
%[text] Both remaining factors are now recovered close to the truth, and both intervals sit well inside their priors instead of running to the bounds.
bar([Tci3.Range(1:2,2)-Tci3.Range(1:2,1), Tci2.Range(:,2)-Tci2.Range(:,1)])
set(gca,'XTickLabel', T2.Properties.RowNames)
ylabel('95% confidence interval width')
legend('all three free','V fixed','Location','northeast')
title('Fixing one factor sharpens the other two')
grid on
%%
%[text] ## 6. What the example shows
%[text] Compare the two runs on the two things people usually conflate:
%[text:table]
%[text] | | all three free | V fixed |
%[text] | --- | --- | --- |
%[text] | cost of the best fit | lower | slightly **higher** |
%[text] | correlation between factors | near $ \\pm 1 $ | much weaker |
%[text] | CI for $ k_a $ | the entire prior | a usable interval |
%[text] | estimates vs truth | biased | close |
%[text:table]
%[text] Fixing $ V $ made the fit *worse* and the science *better*. Cost measures how well a curve passes through points; it does not measure whether the factors that produced that curve could have been recovered. Only the identifiability analysis answers that, which is why it belongs in the workflow rather than after it.
%[text] Where to go next: |gsua_ia| adds correlation heatmaps, fit-quality filtering and detection of multiple global minima across the multistart pool; |gsua_covmetric| scores an uncertainty band for the same confounding signature; and |gsua_dmatrix| with |'Method','Joint'| samples that pool while preserving its correlation structure. The symbolic-math companion to this example works the same problem for an ODE model defined with Symbolic Math Toolbox.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
