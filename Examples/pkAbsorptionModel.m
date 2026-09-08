function sol = pkAbsorptionModel(params, domain, ~)
%PKABSORPTIONMODEL One-compartment pharmacokinetic model with first-order absorption.
%
%   sol = pkAbsorptionModel(params, domain, opt)
%
%   Plasma concentration after a single oral dose D:
%
%       c(t) = D*ka / (V*(ka-ke)) * ( exp(-ke*t) - exp(-ka*t) )
%
%   params = [ka; ke; V]
%     ka <-- absorption rate constant (1/h)
%     ke <-- elimination rate constant (1/h)
%     V  <-- apparent volume of distribution (L)
%   domain = [t0 tf] time span in hours.
%
%   Returns an ODE-solver-shaped struct (sol.x, sol.y) so that GSUA_EVAL can
%   interpolate the solution onto any xdata, exactly as it does for a model
%   integrated with ode45 or dde23. The model is closed-form, so the dense
%   output costs nothing.
%
%   ka is deliberately bounded above ke in the examples that use this model:
%   at ka == ke the expression above is singular (its limit is the finite
%   D*ka*t*exp(-ka*t)/V, but the formula itself is 0/0), and swapping the two
%   leaves c(t) unchanged -- the classic flip-flop ambiguity, which would make
%   the parameters unidentifiable for an uninteresting algebraic reason rather
%   than an experimental-design one.
%
%   The unused third input is the 'opt' configuration slot documented in the
%   user guide's "Time-dependent functions" section. It is required, not
%   optional: gsua_eval calls a time-dependent user function as
%   fun(pars,domain,opt), so a two-input function is rejected at evaluation
%   time even though gsua_userdefined will happily build a handle for one.
%   The toolbox's own Examples/user_dependent.m takes the same three inputs.
%
%   See also GSUA_DATAPREP, GSUA_EVAL.
D = 100;                                     % administered dose (mg), known
ka = params(1);
ke = params(2);
V  = params(3);
sol.x = linspace(domain(1), domain(2), 400);
sol.y = (D*ka)./(V*(ka-ke)) .* (exp(-ke*sol.x) - exp(-ka*sol.x));
end
