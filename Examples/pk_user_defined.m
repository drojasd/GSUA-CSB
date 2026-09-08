%[text] # End-to-end system identification: a user-defined model
%[text] This example takes a pharmacokinetic model from a plain MATLAB function all the way to confidence intervals on its parameters, using a **user-defined** model rather than Simulink or Symbolic Math.
%[text] The point is not that the fit succeeds. It is that **a good fit tells you almost nothing about whether your parameters are identifiable** — and that the toolbox will tell you the difference if you ask it.
%[text] A Python notebook version of this example runs the same computation with the same section numbering, so the two can be read side by side.
%[text] ## The model
%[text] A one-compartment model with first-order absorption, the standard description of an orally administered drug:
%[text] $ c(t)=\\frac{D\\,k_a}{V(k_a-k_e)}\\left(e^{-k_e t}-e^{-k_a t}\\right) $
%[text] Three factors are to be identified, with the dose $ D=100 $ mg known:
%[text] - $ k_a $ — absorption rate constant (1/h)
%[text] - $ k_e $ — elimination rate constant (1/h)
%[text] - $ V $ — apparent volume of distribution (L) \
%[text] The implementation is an ordinary function file. It returns an ODE-solver-shaped struct (`sol.x`, `sol.y`) so that `gsua_eval` can interpolate it onto any set of sampling times, exactly as it would for a model integrated with `ode45`:
type pkAbsorptionModel.m
%%
%[text] ## 1. Preparing the environment
%[text] `gsua_dataprep` builds the summary table `T` that every other toolbox function consumes. For a user-defined model it needs the function name, the factor bounds, and the time domain.
% Factor bounds: one row per factor, [lower upper].
ranges = [0.6  3.0      % ka  absorption rate (1/h)
          0.05 0.5      % ke  elimination rate (1/h)
          5    40];     % V   volume of distribution (L)
% 'domain' is the time span the model is integrated over; 'names'/'out_names' are
% cosmetic but propagate into every plot and table the toolbox produces.
[T,~] = gsua_dataprep('pkAbsorptionModel', ranges, 'domain', [0 24], ...
    'names', {'ka','ke','V'}, 'out_names', {'concentration'});
T
%[text] The $ k_a $ bounds are deliberately kept above the $ k_e $ bounds. At $ k_a=k_e $ the closed form is singular, and swapping the two leaves $ c(t) $ unchanged — the classic *flip-flop* ambiguity. Excluding it keeps this example about experimental design rather than an algebraic accident.
%%
%[text] ## 2. Synthetic data
%[text] Working from synthetic data means the truth is known, so the confidence intervals can be *checked* rather than merely reported.
rng(0,'twister')                          % fix the noise draw so the page reproduces
truth = [1.2; 0.25; 15];                  % the values we will try to recover
T.Nominal = truth;
xdata = [0.25 0.5 1 1.5 2 3 4 6 8 10 12 16 20 24];   % sampling schedule (hours)
% gsua_eval(values, T, xdata, ydata, parallel, show, verbose): the trailing two false
% flags suppress the automatic diagnostic plot and the "Progress: N %" print.
clean = gsua_eval(truth, T, xdata, [], false, false, false);
ydata = clean .* (1 + 0.08*randn(size(clean)));      % 8% proportional measurement noise
tdense = linspace(0.05, 24, 300);                    % dense grid, for drawing curves only
plot(tdense, gsua_eval(truth, T, tdense, [], false, false, false), 'LineWidth', 1.5)
hold on
plot(xdata, ydata, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
hold off
xlabel('time (h)')
ylabel('concentration (mg/L)')
legend('true model','measurements','Location','northeast')
title('Simulated single-dose concentration data')
grid on
%%
%[text] ## 3. Can these parameters be estimated at all?
%[text] Before spending any optimizer budget it is worth asking whether the data is even *reachable*: does it fall inside the range of behaviours the model can produce over the factor bounds declared in section 1? If it does not, no amount of optimization will help — the model structure or the bounds are wrong, and that has to be fixed first. This is the reachability check that opens the toolbox's semi-automated identification cycle.
%[text] `gsua_dmatrix` samples the factor box and `gsua_ua` runs the Monte-Carlo ensemble over those samples. `gsua_ua` also applies Monte-Carlo filtering automatically and plots it, which is what the two figures below show: for each factor, how the low and high halves of its sampled range map onto model output.
M0 = gsua_dmatrix(T, 500);                           % 500 samples of the factor box
Y0 = gsua_ua(M0, T, 'xdata', xdata, 'ynom', ydata, 'parallel', false, 'verbose', false);
%[text] `gsua_covmetric` reduces that ensemble to a 5–95% band and scores it. The containment fraction — how much of the measured data actually falls inside the band — is the number to read at this stage.
[cost_data, cost_band, P5, ~, P95] = gsua_covmetric(Y0, ydata, 'margin', 0.1);
reachable = mean(ydata >= P5 & ydata <= P95);
table(reachable, median(P95-P5), cost_data, cost_band, ...
    'VariableNames', {'contained','median_band_width','cost_data','cost_band'})
%[text] All of the data lies inside the reachable band, so estimation is worth attempting. Note that `cost_data` and `cost_band` are both far above 1 here, and that is expected rather than alarming: they are normalized against a tight tolerance and are meaningful as a *post-convergence* check (section 7), not as a pass/fail gate against a prior range this wide. The band is enormous — a median width of about 5 mg/L against data that never exceeds 4.5 — which is exactly what an uninformative prior looks like before any fitting.
plot(xdata, P5, 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
hold on
plot(xdata, P95, 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
plot(xdata, ydata, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
hold off
xlabel('time (h)')
ylabel('concentration (mg/L)')
legend('5th percentile','95th percentile','measurements','Location','northeast')
title(sprintf('Reachable band before fitting (%.0f%% of data contained)', 100*reachable))
grid on
%%
%[text] ## 4. Estimating the parameters
%[text] `gsua_pe` runs a multistart estimation: `'N',20` restarts the optimizer from twenty different points in the factor space, which is how you find out whether the problem has one solution or several.
% 'margin',0.1 selects the correlation-penalized cost and records the margin on the
% output table, so functions further down the pipeline can recover what was scored.
% 'timer',false suppresses the progress readout.
[T3,res3] = gsua_pe(T, xdata, ydata, 'solver','lsqc', 'N',20, 'margin',0.1, 'timer',false);
Est3 = T3.Estlsqc;                        % 3 x 20: one column per multistart run
table(truth, Est3(:,1), 'VariableNames', {'true','estimated'}, 'RowNames', T3.Properties.RowNames)
%[text] Every one of the twenty runs converged to the same cost, and the fitted curve passes cleanly through the data. On most projects this is where the analysis would stop.
plot(tdense, gsua_eval(Est3(:,1), T3, tdense, [], false, false, false), 'LineWidth', 1.5)
hold on
plot(xdata, ydata, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
hold off
xlabel('time (h)')
ylabel('concentration (mg/L)')
legend('fitted model','measurements','Location','northeast')
title(sprintf('Fit with all three factors free (cost = %.4g)', min(res3)))
grid on
%%
%[text] ## 5. Diagnosing identifiability
%[text] The correlation between the repeated estimates is the first warning sign. Values near $ \\pm 1 $ mean the factors trade off against each other: many different combinations reproduce the same curve.
array2table(corr(Est3'), 'VariableNames', T3.Properties.RowNames, 'RowNames', T3.Properties.RowNames)
%[text] `gsua_likelihood` turns that into something actionable. It profiles each factor — stepping it away from the estimate while re-fitting all the others at every step — and reports the interval over which the fit stays statistically acceptable.
% Positional arguments: (T, xdata, ydata, alpha, step, margin, tol1, tol2, limit,
% reps, show, parallel, saver, pars). alpha = 0.95 is the confidence level; limit is
% the bisection budget per bound; pars = [] profiles every factor.
%
% margin is a RELATIVE STANDARD DEVIATION OFFSET BY ONE here: gsua_likelihood uses
% (margin-1) as the assumed relative noise, so 1.08 means the 8% noise the data
% actually carries. Passing 0.1 would assert 90% noise -- and would also make
% gsua_pe's internal +1 offset positive, silently switching its inner refit from the
% likelihood to plain least squares. The Python port drops the offset, so its
% margin=0.08 is this margin=1.08.
Tci3 = gsua_likelihood(T3, xdata, ydata, 0.95, 0.05, 1.08, 0.01, 0.01, 15, 1, false, false, false, [], false);
table(Tci3.Range(:,1), Tci3.Range(:,2), Tci3.Range(:,2)-Tci3.Range(:,1), ranges(:,1), ranges(:,2), ...
    'VariableNames', {'CI_low','CI_high','width','prior_low','prior_high'}, ...
    'RowNames', T3.Properties.RowNames)
%[text] This is the result worth stopping on. The correlation between $ k_a $ and $ k_e $ is $ -0.998 $: the two rate constants are very nearly a single degree of freedom, and the intervals are correspondingly loose for a fit this good. A model that reproduces the data this well is still telling us that these three factors are barely separable from one oral concentration curve, which is a textbook pharmacokinetic result rather than a failure of the optimizer.
%%
%[text] ## 6. The remedy: fix what another experiment already knows
%[text] The standard resolution is to measure $ V $ separately, in an intravenous study where it *is* directly identifiable, and then estimate only the two rate constants. In this toolbox a factor is fixed by giving it a degenerate range; fixed factors drop out of `T` entirely, so the table returned below has two rows rather than three.
rangesFixed = [0.6 3.0; 0.05 0.5; 15 15];            % V pinned at its known value
[Tf,~] = gsua_dataprep('pkAbsorptionModel', rangesFixed, 'domain', [0 24], ...
    'names', {'ka','ke','V'}, 'out_names', {'concentration'});
Tf.Nominal = truth(1:height(Tf));                    % only the free factors remain
[T2,res2] = gsua_pe(Tf, xdata, ydata, 'solver','lsqc', 'N',20, 'margin',0.1, 'timer',false);
Tci2 = gsua_likelihood(T2, xdata, ydata, 0.95, 0.05, 1.08, 0.01, 0.01, 15, 1, false, false, false, [], false);
table(truth(1:2), T2.Estlsqc(:,1), Tci2.Range(:,1), Tci2.Range(:,2), Tci2.Range(:,2)-Tci2.Range(:,1), ...
    'VariableNames', {'true','estimated','CI_low','CI_high','width'}, ...
    'RowNames', T2.Properties.RowNames)
%[text] Both remaining factors are now recovered close to the truth, and both intervals sit well inside their bounds instead of running to them.
%%
%[text] ## 7. The two runs side by side
%[text] Comparing the two fits on the quantities people usually conflate:
summary = table([min(res3); corr2free(Est3); Tci3.Range(1,2)-Tci3.Range(1,1); Tci3.Range(2,2)-Tci3.Range(2,1)], ...
                [min(res2); corr2free(T2.Estlsqc); Tci2.Range(1,2)-Tci2.Range(1,1); Tci2.Range(2,2)-Tci2.Range(2,1)], ...
    'VariableNames', {'all_three_free','V_fixed'}, ...
    'RowNames', {'best cost','corr(ka,ke)','CI width ka','CI width ke'})
bar([Tci3.Range(1:2,2)-Tci3.Range(1:2,1), Tci2.Range(:,2)-Tci2.Range(:,1)])
set(gca,'XTickLabel', T2.Properties.RowNames)
ylabel('95% confidence interval width')
legend('all three free','V fixed','Location','northeast')
title('Fixing one factor sharpens the other two')
grid on
%%
%[text] ## 8. What this example shows
%[text] Fixing $ V $ made the fit slightly **worse** and the science considerably **better**: the correlation between $ k_a $ and $ k_e $ collapses, both intervals tighten, and the estimates move onto the truth.
%[text] Cost measures how well a curve passes through points. It does not measure whether the factors that produced that curve could have been recovered. Only the identifiability analysis answers that, which is why it belongs *inside* the workflow rather than after it.
%[text] The companion symbolic-math example reaches the same conclusion from the opposite direction: there, the dataset that fits *better* is the one whose parameters are *less* identifiable.
%[text] Where a multistart run does spread across the factor space instead of converging to a single point, `gsua_ia` adds correlation heatmaps, fit-quality filtering and detection of multiple global minima; `gsua_covmetric` scores an uncertainty band for the same confounding signature; and `gsua_dmatrix` with `'Method','Joint'` propagates the accepted set without destroying its correlation structure.
function r = corr2free(E)
% Correlation between the first two factors across the multistart pool.
R = corr(E');
r = R(1,2);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
