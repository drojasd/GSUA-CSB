"""Model abstraction: the Python replacement for the MATLAB toolbox's table+CustomProperties hack.

The MATLAB toolbox stores a model's callable ("Solver"), its domain, its output selection, and
whether each factor is fixed as ``CustomProperties`` bolted onto a plain ``table`` via ``addprop``,
then dispatches on an integer ``Kind`` (1=Simulink, 2=vectorized user function, 3/4=symbolic ODE via
ode45/ode4, 5=user function with domain+opt, 6=scalar user function) inside every downstream
function (``gsua_deval``, ``gsua_pardeval``, ``gsua_eval``, ...). That was a workaround for MATLAB
not having a clean way to attach a callable + metadata to tabular data -- Python has no reason to
inherit it.

Here a :class:`Model` is a small class carrying parameter metadata (names, range, nominal, fixed)
plus an ``evaluate``/``evaluate_batch`` pair. Concrete subclasses (``UserFunctionModel`` here;
``SymbolicODEModel`` in a later module) replace the ``Kind`` dispatch with real polymorphism.
Simulink-backed models (MATLAB ``Kind == 1``) have no Python equivalent and are intentionally not
represented here -- a Simulink model is accessed from Python via the MATLAB Engine API, not
reimplemented.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable, Sequence
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

RangeMethod = Literal["range", "percent", "normal", "std"]


def build_range(values: ArrayLike, spread: ArrayLike | None, method: RangeMethod) -> NDArray[np.float64]:
    """Build an (Np, 2) [lower, upper] range array, mirroring the MATLAB toolbox's ``rMethod``.

    Args:
        values: For ``method="range"``, an (Np, 2) array of [lower, upper] bounds directly (``spread``
            is ignored). For every other method, an (Np,) array of center/nominal values.
        spread: (Np,) array interpreted per ``method``: percent-of-center for ``"percent"``, standard
            deviation for ``"normal"``/``"std"``. Ignored for ``method="range"``.
        method: ``"range"`` -- ``values`` are already [lower, upper] bounds.
            ``"percent"`` -- ``[v*(1-s/100), v*(1+s/100)]``.
            ``"normal"``/``"std"`` -- ``[v-s, v+s]`` (for ``"normal"``, callers typically sample from
            a normal distribution over this interval rather than treating it as a hard bound; see
            :mod:`gsua_csb._sampling`).

    Returns:
        (Np, 2) float64 array of [lower, upper] bounds. ``lower == upper`` marks a fixed factor,
        the same convention ``gsua_dmatrix``/``gsua_sa``/... use throughout the MATLAB toolbox.

    Raises:
        ValueError: If ``method`` is not one of the four supported values, or ``spread`` is missing
            when required.
    """
    values = np.asarray(values, dtype=np.float64)

    if method == "range":
        arr = values
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("method='range' requires an (Np, 2) [lower, upper] array")
        return arr

    if spread is None:
        raise ValueError(f"method={method!r} requires `spread`")
    spread = np.asarray(spread, dtype=np.float64)

    if method == "percent":
        return np.column_stack([values * (1 - spread / 100), values * (1 + spread / 100)])
    if method in ("normal", "std"):
        return np.column_stack([values - spread, values + spread])
    raise ValueError(f"Unknown range method: {method!r}")


class Model(ABC):
    """Base class for a GSUA model: parameter metadata plus how to evaluate it.

    Attributes:
        names: Parameter names, length Np. Row names in the MATLAB table.
        range: (Np, 2) [lower, upper] bounds. ``T.Range`` in MATLAB.
        nominal: (Np,) nominal/reference parameter values. ``T.Nominal`` in MATLAB.
        output_names: Names of the model's output signal(s). ``CustomProperties.Vars`` in MATLAB.
        domain: Evaluation domain (e.g. a time span), or ``None`` for domain-less (scalar) models.
            ``CustomProperties.Domain`` in MATLAB.
        log_scale: (Np,) bool mask -- True where a free parameter should be *sampled and searched*
            in log10 space rather than linear space (no MATLAB equivalent; new capability). Only
            meaningful for free parameters (``~fixed``) with a strictly positive lower bound; a
            fixed parameter's scale is irrelevant since it is never sampled or optimized.
            :func:`gsua_csb.design_matrix` and :func:`gsua_csb.parameter_estimation` both respect
            this. All ``False`` (the previous, only behavior) unless a concrete class sets it
            otherwise -- see :func:`gsua_csb.load_petab`, which populates it automatically from
            PEtab's ``parameterScale`` column, the format's own signal that a parameter spanning
            many orders of magnitude (extremely common for epidemiology/systems-biology rate
            constants) needs log-space sampling to be searchable at all: linear-uniform sampling
            over e.g. ``[1e-13, 1000]`` essentially never lands anywhere near a true value like
            ``2e-12``, and a linear-space local optimizer cannot traverse that distance either.
    """

    names: list[str]
    range: NDArray[np.float64]
    nominal: NDArray[np.float64]
    output_names: list[str]
    domain: NDArray[np.float64] | None
    log_scale: NDArray[np.bool_]

    @property
    def n_params(self) -> int:
        """Number of parameters (Np), including fixed ones."""
        return len(self.names)

    @property
    def fixed(self) -> NDArray[np.bool_]:
        """(Np,) bool mask: True where ``range[:,0] == range[:,1]`` (a fixed factor).

        Equivalent to ``CustomProperties.Fixed`` in MATLAB, computed live from ``range`` here
        instead of stored separately -- there is only ever one source of truth for "is this
        parameter fixed", unlike the MATLAB version where ``gsua_dmatrix`` recomputes and
        re-stores ``Fixed`` as a side effect.
        """
        return self.range[:, 0] == self.range[:, 1]

    def fix(self, name_or_index: str | int, value: float | None = None) -> None:
        """Collapse a factor's range to a point, marking it fixed (in place).

        Mirrors the MATLAB pattern ``T.Range(i,:) = T.Nominal(i)`` used throughout the toolbox
        (e.g. in the semi-automation identification routine's Phase 5) to give up on estimating a
        parameter. Every method here that samples/evaluates already respects ``fixed`` via
        ``range``, so no other bookkeeping is needed.

        Args:
            name_or_index: Parameter name or integer index.
            value: Value to fix at. Defaults to the parameter's current ``nominal`` value.
        """
        i = self.names.index(name_or_index) if isinstance(name_or_index, str) else name_or_index
        v = self.nominal[i] if value is None else value
        self.range[i, :] = v
        self.nominal[i] = v

    @abstractmethod
    def evaluate(self, params: NDArray[np.float64], xdata: ArrayLike | None = None) -> NDArray[np.float64]:
        """Evaluate the model for a single parameter vector.

        Args:
            params: (Np,) parameter values (all Np, including fixed ones -- unlike the MATLAB
                toolbox's ``fixing`` closures, no separate "free-parameters-only" vector is needed
                since numpy indexing makes it just as easy to pass the full vector).
            xdata: Points to evaluate/interpolate at. If ``None``, uses ``self.domain``.

        Returns:
            Model output. Shape depends on the model: (len(xdata),) for a single-output
            time-series model, (len(xdata), n_outputs) for multi-output, or ``(1,)``/scalar for a
            domain-less model.
        """

    def evaluate_batch(
        self, params: NDArray[np.float64], xdata: ArrayLike | None = None
    ) -> NDArray[np.float64]:
        """Evaluate the model for a batch of parameter sets.

        Default implementation loops calling :meth:`evaluate` once per row -- correct for any
        model, but not fast. This is the direct replacement for MATLAB's ``parfor``-based
        ``gsua_pardeval``/``gsua_eval``: subclasses whose underlying callable can accept a batch of
        parameter sets at once (the MATLAB toolbox's ``vectorized=True`` case) should override this
        for real vectorized speed instead of a Python-level loop, which is usually a much bigger win
        than parallelizing the loop itself (no process/thread overhead, no GIL contention).

        Args:
            params: (N, Np) parameter sets, one row per run.
            xdata: Points to evaluate/interpolate at. If ``None``, uses ``self.domain``.

        Returns:
            (N, ...) stacked output, one row per parameter set.
        """
        rows = [self.evaluate(params[i], xdata) for i in range(params.shape[0])]
        return np.stack(rows, axis=0)


class UserFunctionModel(Model):
    """A model wrapping an arbitrary Python callable.

    Replaces MATLAB's ``gsua_userdefined`` (and the ``Kind`` in {2, 5, 6} cases it produces).
    Where MATLAB inspected ``nargin(user_func)`` to decide whether to call ``func(pars)``,
    ``func(pars, domain)``, or ``func(pars, domain, opt)``, this class instead always calls
    ``func(params, domain=self.domain, **opt)`` when ``domain`` is not ``None``, or
    ``func(params)`` when it is -- Python's keyword arguments make the arity-sniffing unnecessary.

    Args:
        func: The model function. Called as ``func(params)`` if ``domain`` is ``None``, otherwise
            ``func(params, domain, **opt)``.
        names: Parameter names, length Np.
        range: (Np, 2) [lower, upper] bounds (build with :func:`build_range` for the
            percent/normal/std input conventions).
        nominal: (Np,) nominal values. Defaults to the midpoint of ``range``.
        domain: Evaluation domain passed to ``func``, or ``None`` for a domain-less (scalar-output)
            function.
        output_names: Names of the output signal(s). Defaults to ``["out"]``.
        vectorized: If True, ``func`` accepts a batch of parameter sets at once (shape (N, Np)) and
            :meth:`evaluate_batch` calls it directly instead of looping. Mirrors MATLAB's
            ``vectorized`` flag (``Kind == 2`` when true).
        opt: Extra keyword arguments passed to ``func`` on every call (mirrors ``CustomProperties.copt``).
        log_scale: (Np,) bool mask, or ``None`` (default: all ``False``) -- see :attr:`Model.log_scale`.
    """

    def __init__(
        self,
        func: Callable[..., NDArray[np.float64]],
        names: Sequence[str],
        range: NDArray[np.float64],
        nominal: NDArray[np.float64] | None = None,
        domain: ArrayLike | None = None,
        output_names: Sequence[str] = ("out",),
        vectorized: bool = False,
        opt: dict | None = None,
        log_scale: ArrayLike | None = None,
    ) -> None:
        self._func = func
        self.names = list(names)
        self.range = np.asarray(range, dtype=np.float64)
        if self.range.shape != (len(self.names), 2):
            raise ValueError(f"range must have shape ({len(self.names)}, 2), got {self.range.shape}")
        self.nominal = (
            np.asarray(nominal, dtype=np.float64)
            if nominal is not None
            else self.range.mean(axis=1)
        )
        self.domain = None if domain is None else np.asarray(domain, dtype=np.float64)
        self.output_names = list(output_names)
        self.vectorized = vectorized
        self.opt = opt or {}
        self.log_scale = (
            np.zeros(len(self.names), dtype=bool)
            if log_scale is None
            else np.asarray(log_scale, dtype=bool)
        )
        if self.log_scale.shape != (len(self.names),):
            raise ValueError(f"log_scale must have shape ({len(self.names)},), got {self.log_scale.shape}")

    def evaluate(self, params: NDArray[np.float64], xdata: ArrayLike | None = None) -> NDArray[np.float64]:
        d = self.domain if xdata is None else np.asarray(xdata, dtype=np.float64)
        if d is None:
            return np.asarray(self._func(params, **self.opt))
        return np.asarray(self._func(params, d, **self.opt))

    def evaluate_batch(
        self, params: NDArray[np.float64], xdata: ArrayLike | None = None
    ) -> NDArray[np.float64]:
        if not self.vectorized:
            return super().evaluate_batch(params, xdata)
        d = self.domain if xdata is None else np.asarray(xdata, dtype=np.float64)
        if d is None:
            return np.asarray(self._func(params, **self.opt))
        return np.asarray(self._func(params, d, **self.opt))

    @classmethod
    def from_bounds(
        cls,
        func: Callable[..., NDArray[np.float64]],
        values: ArrayLike,
        spread: ArrayLike | None = None,
        method: RangeMethod = "range",
        names: Sequence[str] | None = None,
        **kwargs,
    ) -> UserFunctionModel:
        """Convenience constructor building ``range`` via :func:`build_range`.

        Mirrors calling ``gsua_userdefined(func, Range, 'rMethod', method, ...)`` in MATLAB, with
        ``values``/``spread`` playing the role of ``Range`` interpreted according to ``method``
        (see :func:`build_range`).
        """
        rng = build_range(values, spread, method)
        if names is None:
            names = [str(i) for i in range(rng.shape[0])]
        return cls(func, names=names, range=rng, **kwargs)
