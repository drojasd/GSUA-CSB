%[text] # End-to-end system identification: a Symbolic Math model
%[text] This example identifies the parameters of an epidemic model written as a system of differential equations with Symbolic Math Toolbox, and asks the same question its user-defined companion asks: once the model fits, can the parameters actually be recovered?
%[text] Here the answer depends on **when you stopped collecting data**, which makes it a concrete demonstration that identifiability is a property of the experiment and not only of the model.
%[text] A Python notebook version of this example runs the same computation with the same section numbering, so the two can be read side by side.
%[text] ## The model
%[text] The classical SIR system for a closed population of $ N=1000 $:
%[text] $ \\frac{dS}{dt}=-\\beta\\frac{SI}{N},\\qquad \\frac{dI}{dt}=\\beta\\frac{SI}{N}-\\gamma I,\\qquad \\frac{dR}{dt}=\\gamma I $
%[text] Two factors are to be identified:
%[text] - $ \\beta $ — transmission rate (1/day)
%[text] - $ \\gamma $ — recovery rate (1/day) \
%[text] Their ratio is the basic reproduction number $ R_0=\\beta/\\gamma $, the quantity an epidemiologist actually wants. Only the infected compartment $ I(t) $ is observed, and the initial conditions are known: one infected individual in an otherwise susceptible population.
%%
%[text] ## 1. Preparing the environment
%[text] The system is declared with `syms` and handed straight to `gsua_dataprep`, which generates the simulation code and returns the summary table.
syms S(t) I(t) R(t) beta gamma
N = 1000;                                  % closed population
odes = [diff(S) == -beta*S*I/N             % susceptible
        diff(I) ==  beta*S*I/N - gamma*I   % infected  (this is what we observe)
        diff(R) ==  gamma*I];              % recovered
vars = [S I R];                            % state order -- fixes the factor order below
%[text] For a symbolic model the factor vector is the **initial conditions of the states first**, in the order given by `vars`, followed by the model parameters. Giving a factor a degenerate range fixes it and drops it from the table, so the five entries below leave only $ \\beta $ and $ \\gamma $ free.
% Rows: S0, I0, R0 (all known -> degenerate), then beta, gamma (to estimate).
% 'output',2 declares that the SECOND state, I(t), is the measured one.
[T,~] = gsua_dataprep(odes, vars, [0 80], 'sirEpidemicModel', ...
    'range', [999 999      % S0  susceptible at t=0
              1   1        % I0  infected at t=0
              0   0        % R0  recovered at t=0
              0.15 0.9     % beta
              0.02 0.4], 'output', 2);
T
%%
%[text] ## 2. Synthetic data
%[text] The true epidemic uses $ \\beta=0.35 $ and $ \\gamma=0.10 $, so $ R_0=3.5 $. Case counts are noisier when they are larger, so the measurements carry Poisson-like noise whose spread grows with the count.
rng(3,'twister')                           % fix the noise draw so the page reproduces
truth = [0.35; 0.10];                      % beta, gamma -- the values to recover
T.Nominal = truth;
xfull = linspace(1, 80, 20);               % 80 days of surveillance, 20 reports
% Trailing false,false suppresses gsua_eval's automatic diagnostic plot.
cleanFull = gsua_eval(truth, T, xfull, [], false, false);
yfull = cleanFull + sqrt(max(cleanFull,1)).*randn(size(cleanFull));
tdense = linspace(1, 80, 300);             % dense grid, for drawing curves only
plot(tdense, gsua_eval(truth, T, tdense, [], false, false), 'LineWidth', 1.5)
hold on
plot(xfull, yfull, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
xline(25, '--', 'end of early phase', 'LabelVerticalAlignment', 'bottom')
hold off
xlabel('time (days)')
ylabel('infected individuals I(t)')
legend('true epidemic','surveillance data','Location','northeast')
title('Simulated epidemic, R_0 = 3.5')
grid on
%%
%[text] ## 3. Can these parameters be estimated at all?
%[text] Before spending any optimizer budget it is worth asking whether the data is even *reachable*: does it fall inside the range of epidemics the model can produce over the $ \\beta $ and $ \\gamma $ bounds declared in section 1? If it does not, no amount of optimization will help — the model structure or the bounds are wrong, and that has to be fixed first. This is the reachability check that opens the toolbox's semi-automated identification cycle.
%[text] `gsua_dmatrix` samples the factor box and `gsua_ua` runs the Monte-Carlo ensemble over those samples, applying Monte-Carlo filtering and plotting it automatically — that is what the figures below show, factor by factor.
M0 = gsua_dmatrix(T, 300);                          % 300 samples of the (beta, gamma) box
Y0 = gsua_ua(M0, T, 'xdata', xfull, 'ynom', yfull, 'parallel', false);
%[text] `gsua_covmetric` reduces the ensemble to a 5–95% band. The containment fraction is the number to read at this stage; `cost_data` and `cost_band` are normalized against a tight tolerance and only become meaningful after convergence, in section 7.
[cost_data, cost_band, P5, ~, P95] = gsua_covmetric(Y0, yfull, 'margin', 0.1);
reachable = mean(yfull >= P5 & yfull <= P95);
table(reachable, median(P95-P5), cost_data, cost_band, ...
    'VariableNames', {'contained','median_band_width','cost_data','cost_band'})
%[text] Nearly all the observations fall inside the reachable band, so an epidemic of this shape is within the model's declared range and estimation is worth attempting. The band itself is far too wide to be useful as an answer — it spans most of the population — which is precisely the uncertainty the next sections set out to reduce.
plot(xfull, P5, 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
hold on
plot(xfull, P95, 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
plot(xfull, yfull, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
hold off
xlabel('time (days)')
ylabel('infected individuals I(t)')
legend('5th percentile','95th percentile','surveillance data','Location','northeast')
title(sprintf('Reachable band before fitting (%.0f%% of data contained)', 100*reachable))
grid on
%%
%[text] ## 4. Estimating from the full epidemic
%[text] With the whole curve in hand — growth, peak and decline — both factors are estimated by multistart least squares, twenty restarts as before.
[Tfull,resFull] = gsua_pe(T, xfull, yfull, 'solver','lsqc', 'N',20, 'margin',0.1, 'timer',false);
Efull = Tfull.Estlsqc;                     % 2 x 20: one column per multistart run
% margin here is a relative standard deviation OFFSET BY ONE: 1.1 asserts 10% noise.
% Passing 0.1 would assert 90% noise and would also flip gsua_pe's internal +1 offset
% positive, silently switching its inner refit from the likelihood to plain least squares.
TciFull = gsua_likelihood(Tfull, xfull, yfull, 0.95, 0.05, 1.1, 0.01, 0.01, 15, 1, false, false, false, []);
table(truth, Efull(:,1), TciFull.Range(:,1), TciFull.Range(:,2), TciFull.Range(:,2)-TciFull.Range(:,1), ...
    'VariableNames', {'true','estimated','CI_low','CI_high','width'}, ...
    'RowNames', Tfull.Properties.RowNames)
%[text] Both factors land close to the truth and both intervals are narrow. Reading off the quantity that matters:
R0_full = Efull(1,1)/Efull(2,1)
%%
%[text] ## 5. Estimating from the early phase only
%[text] Now suppose the analysis had to be done *during* the outbreak, with only the first 25 days available — the situation every real-time epidemic assessment faces. Nothing about the model changes; only the observation window shrinks.
% Same system, shorter domain. The model file is regenerated for the new time span.
[Te,~] = gsua_dataprep(odes, vars, [0 25], 'sirEpidemicModel', ...
    'range', [999 999; 1 1; 0 0; 0.15 0.9; 0.02 0.4], 'output', 2);
Te.Nominal = truth;
xearly = linspace(1, 25, 12);              % 25 days, 12 reports
cleanEarly = gsua_eval(truth, Te, xearly, [], false, false);
yearly = cleanEarly + sqrt(max(cleanEarly,1)).*randn(size(cleanEarly));
[Tearly,resEarly] = gsua_pe(Te, xearly, yearly, 'solver','lsqc', 'N',20, 'margin',0.1, 'timer',false);
Eearly = Tearly.Estlsqc;
TciEarly = gsua_likelihood(Tearly, xearly, yearly, 0.95, 0.05, 1.1, 0.01, 0.01, 15, 1, false, false, false, []);
table(truth, Eearly(:,1), TciEarly.Range(:,1), TciEarly.Range(:,2), TciEarly.Range(:,2)-TciEarly.Range(:,1), ...
    'VariableNames', {'true','estimated','CI_low','CI_high','width'}, ...
    'RowNames', Tearly.Properties.RowNames)
%[text] The estimates still look plausible. Note the cost, though: the early-phase fit is *better* than the full-epidemic fit was — fewer points, all of them on a smooth exponential rise — while the intervals have widened and $ \\gamma $ now reaches nearly to its lower bound.
%%
%[text] ## 6. Diagnosing identifiability
%[text] The correlation between the repeated estimates explains what happened.
table(corrBetaGamma(Efull), corrBetaGamma(Eearly), ...
    'VariableNames', {'full_epidemic','early_phase'}, 'RowNames', {'corr(beta,gamma)'})
%[text] In the early phase it is essentially $ 1 $. During exponential growth the data constrain only the growth rate, roughly $ \\beta-\\gamma $, so any pair with the right difference reproduces the observations equally well. It takes the peak — where susceptibles are depleted and the curve turns over — to separate the two.
%%
%[text] ## 7. The two windows side by side
summary = table([min(resFull); corrBetaGamma(Efull); TciFull.Range(1,2)-TciFull.Range(1,1); ...
                 TciFull.Range(2,2)-TciFull.Range(2,1); Efull(1,1)/Efull(2,1)], ...
                [min(resEarly); corrBetaGamma(Eearly); TciEarly.Range(1,2)-TciEarly.Range(1,1); ...
                 TciEarly.Range(2,2)-TciEarly.Range(2,1); Eearly(1,1)/Eearly(2,1)], ...
    'VariableNames', {'full_epidemic','early_phase'}, ...
    'RowNames', {'best cost','corr(beta,gamma)','CI width beta','CI width gamma','R0 (true 3.5)'})
plot(tdense, gsua_eval(Efull(:,1), Tfull, tdense, [], false, false), 'LineWidth', 1.5)
hold on
plot(tdense, gsua_eval(Eearly(:,1), Tfull, tdense, [], false, false), '--', 'LineWidth', 1.5)
plot(xfull, yfull, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
hold off
xlabel('time (days)')
ylabel('infected individuals I(t)')
legend('fit to full epidemic','fit to early phase, extrapolated','data','Location','northeast')
title('Where the early-phase fit leads')
grid on
%%
%[text] ## 8. What this example shows
%[text] The same model and the same estimator produced two very different states of knowledge, and the difference was the observation window rather than anything about the algorithm. The better-fitting dataset was the less informative one.
%[text] Fit quality measures agreement with the points you have; identifiability measures whether those points could have distinguished your parameters from the alternatives. They are different questions, and only the second tells you whether an estimate is worth reporting.
%[text] For a real outbreak the consequence follows directly: an $ R_0 $ estimated before the peak carries a confidence interval wide enough to change policy conclusions, and quoting the point estimate alone would hide that.
%[text] The companion user-defined example reaches the same conclusion from the opposite direction: there, adding information by *fixing* a known factor makes the fit slightly worse and the parameters recoverable. Where a multistart run does spread across the factor space, `gsua_ia` adds correlation heatmaps and detection of multiple global minima, `gsua_noisefloor` calibrates which fits to accept against the observation noise, and `gsua_dmatrix` with `'Method','Joint'` propagates the accepted set without destroying its correlation structure.
function r = corrBetaGamma(E)
% Correlation between beta and gamma across the multistart pool.
R = corr(E');
r = R(1,2);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
