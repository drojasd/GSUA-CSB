"""Symbolic-ODE models: SymPy for the system, ``scipy.integrate.solve_ivp`` for numeric solving.

Python port of ``gsua_dpmat`` + ``gsua_odefun``. Replaces MATLAB's mass-matrix extraction and the
bundled fixed-step ``ode4`` Runge-Kutta integrator: ``solve_ivp`` already provides both adaptive
(``"RK45"``, the ``ode45`` analog) and other methods, and natively evaluates at arbitrary
``t_eval`` points -- so there is no separate ``gsua_intrp``-style interpolation step needed here.

Scope: explicit ODEs (``dy/dt = f(t, y, params)``) only -- MATLAB's implicit/DAE support via
``massMatrixForm`` is not ported. ``solve_ivp`` does not solve DAEs directly; supporting them would
need a different integrator entirely and is a separate, larger addition, not a gap in this class.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence

import numpy as np
import sympy as sp
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp

from ._model import Model, RangeMethod, build_range


class SymbolicODEModel(Model):
    """An ODE model defined symbolically, solved numerically via ``scipy.integrate.solve_ivp``.

    Mirrors ``gsua_dpmat``: the factor vector (``params`` everywhere else in this package) is the
    concatenation of ``[initial conditions of the state variables, true model parameters]``, in
    that order -- the same convention the MATLAB toolbox uses so sensitivity analysis, parameter
    estimation, etc. treat initial conditions and parameters uniformly as "factors" to sample or
    estimate.

    Args:
        odes: RHS expressions, one per state variable: ``d(state_vars[i])/dt = odes[i]``.
        state_vars: SymPy symbols for the state variables, in the same order as ``odes``.
        t: The SymPy symbol used for time in ``odes``.
        params: SymPy symbols for the true (non-initial-condition) parameters appearing in
            ``odes``, in the order they should appear in the factor vector after the initial
            conditions.
        domain: Default evaluation domain (time points), used when ``evaluate``/``evaluate_batch``
            are called without an explicit ``xdata``.
        names: Factor names, length ``len(state_vars) + len(params)`` (initial conditions first,
            then parameters). Defaults to the SymPy symbols' names.
        range: (Np, 2) bounds, same convention as :class:`gsua_csb.Model`. Build with
            :func:`gsua_csb.build_range` for the percent/normal/std input conventions, or use
            :meth:`from_bounds`.
        nominal: (Np,) nominal values. Defaults to the midpoint of ``range``.
        method: Passed to ``solve_ivp``. Default ``"RK45"`` (the ``ode45`` analog); use
            ``"Radau"``/``"BDF"`` for stiff systems.
        solver_kwargs: Extra keyword arguments passed to every ``solve_ivp`` call (e.g. ``rtol``,
            ``atol``, ``max_step``).

    Note:
        :meth:`evaluate_batch` uses the base class's default loop -- each run is an independent
        initial value problem, so there is no batch-vectorization equivalent to
        ``UserFunctionModel``'s ``vectorized=True`` path here. For real speedup across many runs,
        parallelize externally (e.g. ``joblib.Parallel``) rather than expecting this class to do it.
    """

    def __init__(
        self,
        odes: Sequence[sp.Expr],
        state_vars: Sequence[sp.Symbol],
        t: sp.Symbol,
        params: Sequence[sp.Symbol],
        domain: ArrayLike,
        names: Sequence[str] | None = None,
        range: NDArray[np.float64] | None = None,
        nominal: NDArray[np.float64] | None = None,
        method: str = "RK45",
        solver_kwargs: dict | None = None,
    ) -> None:
        self.n_states = len(state_vars)
        self.n_true_params = len(params)
        n = self.n_states + self.n_true_params

        self._rhs = sp.lambdify((t, list(state_vars), list(params)), list(odes), modules="numpy")
        self.method = method
        self.solver_kwargs = solver_kwargs or {}

        self.names = (
            list(names)
            if names is not None
            else [str(s) for s in state_vars] + [str(p) for p in params]
        )
        if len(self.names) != n:
            raise ValueError(f"names must have length {n} (n_states + n_params), got {len(self.names)}")

        self.range = np.asarray(range, dtype=np.float64) if range is not None else np.zeros((n, 2))
        if self.range.shape != (n, 2):
            raise ValueError(f"range must have shape ({n}, 2), got {self.range.shape}")
        self.nominal = (
            np.asarray(nominal, dtype=np.float64) if nominal is not None else self.range.mean(axis=1)
        )
        self.domain = np.asarray(domain, dtype=np.float64)
        self.output_names = list(self.names[: self.n_states])

    def evaluate(self, params: NDArray[np.float64], xdata: ArrayLike | None = None) -> NDArray[np.float64]:
        d = self.domain if xdata is None else np.asarray(xdata, dtype=np.float64)
        y0 = params[: self.n_states]
        true_params = params[self.n_states :]

        sol = solve_ivp(
            lambda tt, yy: self._rhs(tt, yy, true_params),
            (d[0], d[-1]),
            y0,
            t_eval=d,
            method=self.method,
            **self.solver_kwargs,
        )
        if not sol.success:
            warnings.warn(f"solve_ivp failed: {sol.message}", stacklevel=2)
            return np.full((self.n_states, len(d)), np.inf)
        return sol.y

    @classmethod
    def from_bounds(
        cls,
        odes: Sequence[sp.Expr],
        state_vars: Sequence[sp.Symbol],
        t: sp.Symbol,
        params: Sequence[sp.Symbol],
        domain: ArrayLike,
        values: ArrayLike,
        spread: ArrayLike | None = None,
        method: RangeMethod = "range",
        names: Sequence[str] | None = None,
        **kwargs,
    ) -> SymbolicODEModel:
        """Convenience constructor building ``range`` via :func:`gsua_csb.build_range`.

        Mirrors calling ``gsua_dpmat(odes, vars, domain, name, 'rMethod', method, ...)`` in MATLAB.
        """
        rng = build_range(values, spread, method)
        return cls(odes, state_vars, t, params, domain, names=names, range=rng, **kwargs)
