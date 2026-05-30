#!/usr/bin/env python3
"""Ajuste riguroso de la característica I-V de un diodo.

El script compara varios modelos físicamente distintos y evita el problema
común de ajustar V(I) con una resistencia en paralelo mal planteada.  Para un
medidor experimental normalmente se controla/lee el voltaje y se mide la
corriente, por eso se ajusta I(V).

Modelos implementados:
  M0: offset instrumental              I = Ioff
  M1: Shockley ideal                   I = Ioff + Is*(exp(V/(nVt))-1)
  M2: Shockley + Rs                    I = Ioff + solución explícita con Lambert W
  M3: Shockley + Rp                    I = Ioff + Is*(exp(V/(nVt))-1) + V/Rp
  M4: Shockley + Rs + Rp               I = Ioff + solución implícita robusta

Uso:
  ./ajuste_diodo_riguroso.py

Edita la sección CONFIGURACIÓN si tu CSV, columnas o unidades tienen otros nombres.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares
from scipy.special import lambertw

warnings.filterwarnings("ignore", category=RuntimeWarning)

K_B = 1.380649e-23
Q_E = 1.602176634e-19
T_REF = 300.0
VT_REF = K_B * T_REF / Q_E

# ========================== CONFIGURACIÓN DEL USUARIO =========================
# Este archivo está pensado como ejecutable directo, no como CLI configurable.
# Cambia estos valores aquí si tu archivo o tus columnas tienen otros nombres.
CSV_PATH = Path("datos.csv")
VOLTAGE_COL = "Voltaje"
CURRENT_COL = "Corriente"
CURRENT_UNIT = "mA"  # opciones: "A", "mA", "uA", "µA", "nA"
PLOT_PATH = Path("ajuste_diodo_riguroso.png")
GENERAR_GRAFICO = True
EXPORT_THEORY_CSV = Path("valores_teoricos_modelos_20uA_100uA.csv")
CORRIENTE_MIN_EXPORT_A = 20e-6
CORRIENTE_MAX_EXPORT_A = 100e-6
PUNTOS_EXPORT = 200
# =============================================================================


@dataclass(frozen=True)
class FitSpec:
    name: str
    param_names: tuple[str, ...]
    p0: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    current_fn: Callable[[np.ndarray, np.ndarray], np.ndarray]


@dataclass
class FitResult:
    spec: FitSpec
    params: np.ndarray
    stderr: np.ndarray
    current: np.ndarray
    residuals: np.ndarray
    rmse_a: float
    rmse_log: float
    r2: float
    aic: float
    bic: float
    success: bool
    message: str


def safe_expm1(x: np.ndarray | float) -> np.ndarray | float:
    """expm1 con recorte para evitar overflow sin ocultar la física del modelo."""
    return np.expm1(np.clip(x, -745.0, 700.0))


def current_offset(v: np.ndarray, p: np.ndarray) -> np.ndarray:
    (ioff,) = p
    return np.full_like(v, ioff, dtype=float)


def current_ideal(v: np.ndarray, p: np.ndarray) -> np.ndarray:
    log_is, n, ioff = p
    isat = np.exp(log_is)
    a = n * VT_REF
    return ioff + isat * safe_expm1(v / a)


def current_rs(v: np.ndarray, p: np.ndarray) -> np.ndarray:
    log_is, n, log_rs, ioff = p
    isat = np.exp(log_is)
    rs = np.exp(log_rs)
    a = n * VT_REF
    z_log = np.log(rs * isat / a) + (v + rs * isat) / a
    z = np.exp(np.clip(z_log, -745.0, 700.0))
    diode_current = (a / rs) * np.real(lambertw(z)) - isat
    return ioff + diode_current


def current_rp(v: np.ndarray, p: np.ndarray) -> np.ndarray:
    log_is, n, log_rp, ioff = p
    isat = np.exp(log_is)
    rp = np.exp(log_rp)
    a = n * VT_REF
    return ioff + isat * safe_expm1(v / a) + v / rp


def _implicit_current_one(v: float, isat: float, a: float, rs: float, rp: float) -> float:
    """Resuelve I = Is*(exp((V-I*Rs)/a)-1) + (V-I*Rs)/Rp."""

    def f(i: float) -> float:
        vd = v - i * rs
        return isat * safe_expm1(vd / a) + vd / rp - i

    # El circuito es monótono para Rs>=0 y Rp>0; expandimos un intervalo amplio.
    scale = max(abs(v) / max(rs, 1e-30), abs(v) / rp, isat, 1e-12)
    lo = -10.0 * scale - 1e-9
    hi = 10.0 * scale + 1e-9
    flo = f(lo)
    fhi = f(hi)
    for _ in range(80):
        if np.isfinite(flo) and np.isfinite(fhi) and flo * fhi <= 0:
            return brentq(f, lo, hi, maxiter=100, xtol=1e-14, rtol=1e-12)
        lo *= 2.0
        hi *= 2.0
        flo = f(lo)
        fhi = f(hi)
    raise RuntimeError(f"No se pudo acotar la raíz implícita para V={v:g}")


def current_rs_rp(v: np.ndarray, p: np.ndarray) -> np.ndarray:
    log_is, n, log_rs, log_rp, ioff = p
    isat = np.exp(log_is)
    rs = np.exp(log_rs)
    rp = np.exp(log_rp)
    a = n * VT_REF
    return ioff + np.array([_implicit_current_one(float(x), isat, a, rs, rp) for x in v])


def asinh_residual(model_i: np.ndarray, data_i: np.ndarray, i_scale: float) -> np.ndarray:
    """Residuo casi relativo para corrientes grandes y lineal cerca de cero."""
    return np.arcsinh(model_i / i_scale) - np.arcsinh(data_i / i_scale)


def finite_metrics(v: np.ndarray, i: np.ndarray, pred: np.ndarray, resid: np.ndarray, k: int) -> tuple[float, float, float, float, float]:
    n_obs = len(i)
    rmse_a = float(np.sqrt(np.mean((pred - i) ** 2)))
    rmse_log = float(np.sqrt(np.mean(resid**2)))
    ss_res = float(np.sum((pred - i) ** 2))
    ss_tot = float(np.sum((i - np.mean(i)) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    rss = max(float(np.sum(resid**2)), np.finfo(float).tiny)
    aic = float(n_obs * np.log(rss / n_obs) + 2 * k)
    bic = float(n_obs * np.log(rss / n_obs) + k * np.log(n_obs))
    return rmse_a, rmse_log, r2, aic, bic


def covariance_from_jacobian(jac: np.ndarray, residuals: np.ndarray) -> np.ndarray:
    dof = max(1, jac.shape[0] - jac.shape[1])
    s_sq = float(np.sum(residuals**2) / dof)
    _, s, vt = np.linalg.svd(jac, full_matrices=False)
    threshold = np.finfo(float).eps * max(jac.shape) * s[0]
    s = s[s > threshold]
    vt = vt[: len(s)]
    return (vt.T / s**2) @ vt * s_sq


def fit_model(spec: FitSpec, v: np.ndarray, i: np.ndarray, i_scale: float) -> FitResult:
    starts = [spec.p0]
    rng = np.random.default_rng(20260530)
    for _ in range(24):
        jitter = rng.normal(0.0, 0.35, size=spec.p0.size)
        starts.append(np.clip(spec.p0 + jitter, spec.lower, spec.upper))

    best = None
    best_cost = np.inf

    def residual(p: np.ndarray) -> np.ndarray:
        pred = spec.current_fn(v, p)
        if not np.all(np.isfinite(pred)):
            return np.full_like(i, 1e12, dtype=float)
        return asinh_residual(pred, i, i_scale)

    for p0 in starts:
        try:
            opt = least_squares(
                residual,
                p0,
                bounds=(spec.lower, spec.upper),
                loss="soft_l1",
                f_scale=1.0,
                x_scale="jac",
                max_nfev=20000,
            )
        except Exception:
            continue
        if opt.cost < best_cost:
            best = opt
            best_cost = opt.cost

    if best is None:
        nan = np.full(spec.p0.size, np.nan)
        return FitResult(spec, nan, nan, np.full_like(i, np.nan), np.full_like(i, np.nan), np.nan, np.nan, np.nan, np.nan, np.nan, False, "falló")

    pred = spec.current_fn(v, best.x)
    resid = asinh_residual(pred, i, i_scale)
    rmse_a, rmse_log, r2, aic, bic = finite_metrics(v, i, pred, resid, len(best.x))
    try:
        cov = covariance_from_jacobian(best.jac, resid)
        stderr = np.sqrt(np.diag(cov))
    except Exception:
        stderr = np.full_like(best.x, np.nan)

    return FitResult(spec, best.x, stderr, pred, resid, rmse_a, rmse_log, r2, aic, bic, best.success, best.message)


def robust_initial_guess(v: np.ndarray, i: np.ndarray) -> tuple[float, float, float, float, float]:
    positive = i > max(np.nanmax(np.abs(i)) * 1e-9, 1e-12)
    if np.count_nonzero(positive) >= 3:
        vp = v[positive]
        ip = i[positive]
        mid = (ip > np.percentile(ip, 20)) & (ip < np.percentile(ip, 80))
        if np.count_nonzero(mid) >= 2:
            slope, intercept = np.polyfit(vp[mid], np.log(ip[mid]), 1)
            n0 = float(np.clip(1.0 / (slope * VT_REF), 0.8, 8.0)) if slope > 0 else 2.0
            is0 = float(np.clip(np.exp(intercept), 1e-18, 1e-3))
        else:
            n0, is0 = 2.0, 1e-12
        high = ip >= np.percentile(ip, 80)
        if np.count_nonzero(high) >= 2:
            rs0 = abs(float(np.polyfit(ip[high], vp[high], 1)[0]))
        else:
            rs0 = 1.0
    else:
        n0, is0, rs0 = 2.0, 1e-12, 1.0

    low = np.abs(v) <= max(0.1, np.percentile(np.abs(v), 25))
    if np.count_nonzero(low) >= 2:
        conductance = abs(float(np.polyfit(v[low], i[low], 1)[0]))
        rp0 = 1.0 / max(conductance, 1e-12)
        ioff0 = float(np.median(i[low] - conductance * v[low]))
    else:
        rp0 = 1e9
        ioff0 = float(np.median(i))

    return is0, n0, max(rs0, 1e-6), max(rp0, 1.0), ioff0


def build_specs(v: np.ndarray, i: np.ndarray) -> list[FitSpec]:
    is0, n0, rs0, rp0, ioff0 = robust_initial_guess(v, i)
    log_is0 = np.log(is0)
    log_rs0 = np.log(rs0)
    log_rp0 = np.log(rp0)
    i_abs = max(float(np.nanmax(np.abs(i))), 1e-12)

    return [
        FitSpec("M0 offset instrumental", ("Ioff",), np.array([ioff0]), np.array([-10 * i_abs]), np.array([10 * i_abs]), current_offset),
        FitSpec("M1 Shockley ideal", ("ln(Is)", "n", "Ioff"), np.array([log_is0, n0, ioff0]), np.array([np.log(1e-30), 0.5, -10 * i_abs]), np.array([np.log(1e-1), 20.0, 10 * i_abs]), current_ideal),
        FitSpec("M2 Shockley + Rs", ("ln(Is)", "n", "ln(Rs)", "Ioff"), np.array([log_is0, n0, log_rs0, ioff0]), np.array([np.log(1e-30), 0.5, np.log(1e-9), -10 * i_abs]), np.array([np.log(1e-1), 20.0, np.log(1e9), 10 * i_abs]), current_rs),
        FitSpec("M3 Shockley + Rp", ("ln(Is)", "n", "ln(Rp)", "Ioff"), np.array([log_is0, n0, log_rp0, ioff0]), np.array([np.log(1e-30), 0.5, np.log(1.0), -10 * i_abs]), np.array([np.log(1e-1), 20.0, np.log(1e15), 10 * i_abs]), current_rp),
        FitSpec("M4 Shockley + Rs + Rp", ("ln(Is)", "n", "ln(Rs)", "ln(Rp)", "Ioff"), np.array([log_is0, n0, log_rs0, log_rp0, ioff0]), np.array([np.log(1e-30), 0.5, np.log(1e-9), np.log(1.0), -10 * i_abs]), np.array([np.log(1e-1), 20.0, np.log(1e9), np.log(1e15), 10 * i_abs]), current_rs_rp),
    ]


def display_param(name: str, value: float, error: float) -> tuple[str, float, float]:
    if name == "ln(Is)":
        return "Is [A]", float(np.exp(value)), float(abs(np.exp(value) * error))
    if name == "ln(Rs)":
        return "Rs [Ω]", float(np.exp(value)), float(abs(np.exp(value) * error))
    if name == "ln(Rp)":
        return "Rp [Ω]", float(np.exp(value)), float(abs(np.exp(value) * error))
    return name, float(value), float(error)


def load_data(path: Path, voltage_col: str, current_col: str, current_unit: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path)
    missing = {voltage_col, current_col} - set(df.columns)
    if missing:
        raise ValueError(f"Columnas faltantes en {path}: {sorted(missing)}. Columnas disponibles: {list(df.columns)}")
    v = pd.to_numeric(df[voltage_col], errors="coerce").to_numpy(float)
    i = pd.to_numeric(df[current_col], errors="coerce").to_numpy(float)
    factors = {"A": 1.0, "mA": 1e-3, "uA": 1e-6, "µA": 1e-6, "nA": 1e-9}
    i *= factors[current_unit]
    mask = np.isfinite(v) & np.isfinite(i)
    v = v[mask]
    i = i[mask]
    order = np.argsort(v)
    return v[order], i[order]


def make_plot(v: np.ndarray, i: np.ndarray, results: list[FitResult], out: Path) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    ax = axes[0, 0]
    ax.plot(v, i * 1e3, "ko", label="Datos", ms=5)
    for r in results:
        if r.success and r.spec.name != "M0 offset instrumental":
            ax.plot(v, r.current * 1e3, lw=2, label=f"{r.spec.name} (BIC={r.bic:.1f})")
    ax.set_xlabel("Voltaje [V]")
    ax.set_ylabel("Corriente [mA]")
    ax.set_title("Ajuste I(V) en escala lineal")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[0, 1]
    sign = np.sign(i)
    positive = i > 0
    ax.semilogy(v[positive], i[positive], "ko", label="Datos positivos", ms=5)
    for r in results:
        if r.success and r.spec.name != "M0 offset instrumental":
            positive_pred = r.current > 0
            ax.semilogy(v[positive_pred], r.current[positive_pred], lw=2, label=r.spec.name.split()[0])
    ax.set_xlabel("Voltaje [V]")
    ax.set_ylabel("Corriente positiva [A]")
    ax.set_title("Región exponencial")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[1, 0]
    for r in results:
        if r.success and r.spec.name != "M0 offset instrumental":
            ax.plot(v, (r.current - i) * 1e6, "o-", ms=4, label=r.spec.name.split()[0])
    ax.axhline(0, color="k", ls="--", lw=1)
    ax.set_xlabel("Voltaje [V]")
    ax.set_ylabel("Residuo [µA]")
    ax.set_title("Residuos lineales")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    names = [r.spec.name.split()[0] for r in results if r.success]
    bic = [r.bic for r in results if r.success]
    ax.bar(names, bic)
    ax.set_ylabel("BIC menor = mejor")
    ax.set_title("Comparación penalizada por complejidad")
    ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")


def model_slug(name: str) -> str:
    return name.split()[0].lower()


def voltage_for_target_current(result: FitResult, target_i: float, v_data: np.ndarray) -> float:
    """Invierte numéricamente I(V) para exportar V teórico a una corriente dada."""
    if not result.success or result.spec.name == "M0 offset instrumental":
        return np.nan

    def objective(voltage: float) -> float:
        return float(result.spec.current_fn(np.array([voltage]), result.params)[0] - target_i)

    span = max(float(np.ptp(v_data)), 1.0)
    lo = float(np.min(v_data) - 0.25 * span)
    hi = float(np.max(v_data) + 0.25 * span)
    flo = objective(lo)
    fhi = objective(hi)

    for _ in range(60):
        if np.isfinite(flo) and np.isfinite(fhi) and flo * fhi <= 0:
            return float(brentq(objective, lo, hi, maxiter=100, xtol=1e-12, rtol=1e-10))
        lo -= span
        hi += span
        span *= 1.5
        flo = objective(lo)
        fhi = objective(hi)

    return np.nan


def export_theoretical_values(v: np.ndarray, results: list[FitResult], out: Path) -> None:
    """Exporta voltajes teóricos de todos los modelos entre 20 µA y 100 µA."""
    i_min = min(CORRIENTE_MIN_EXPORT_A, CORRIENTE_MAX_EXPORT_A)
    i_max = max(CORRIENTE_MIN_EXPORT_A, CORRIENTE_MAX_EXPORT_A)
    currents = np.linspace(i_min, i_max, PUNTOS_EXPORT)
    table: dict[str, np.ndarray] = {
        "Corriente_objetivo_A": currents,
        "Corriente_objetivo_uA": currents * 1e6,
    }

    for result in sorted(results, key=lambda r: r.spec.name):
        slug = model_slug(result.spec.name)
        table[f"{slug}_Voltaje_teorico_V"] = np.array([
            voltage_for_target_current(result, float(current), v) for current in currents
        ])

    pd.DataFrame(table).to_csv(out, index=False)


def main() -> None:
    v, i = load_data(CSV_PATH, VOLTAGE_COL, CURRENT_COL, CURRENT_UNIT)
    if len(v) < 6:
        raise ValueError("Se necesitan al menos 6 puntos válidos para comparar modelos.")

    i_scale = max(np.percentile(np.abs(i), 10), np.nanmax(np.abs(i)) * 1e-9, 1e-12)
    print("=" * 78)
    print("DIAGNÓSTICO")
    print("=" * 78)
    print(f"Archivo: {CSV_PATH}")
    print(f"Columnas: V='{VOLTAGE_COL}', I='{CURRENT_COL}' ({CURRENT_UNIT})")
    print(f"Puntos válidos: {len(v)}")
    print(f"V: {v.min():.6g} a {v.max():.6g} V")
    print(f"I: {i.min():.6g} a {i.max():.6g} A")
    print(f"Residuo usado: asinh(I/{i_scale:.3e}), estable cerca de cero y por décadas")
    print(f"Vt de referencia: {VT_REF:.6g} V a {T_REF:.1f} K; el parámetro libre es n")

    results = [fit_model(spec, v, i, i_scale) for spec in build_specs(v, i)]
    results.sort(key=lambda r: r.bic if np.isfinite(r.bic) else np.inf)

    print("\n" + "=" * 78)
    print("RESULTADOS ORDENADOS POR BIC")
    print("=" * 78)
    for r in results:
        status = "OK" if r.success else "FALLO"
        print(f"\n{r.spec.name} [{status}]")
        print(f"  RMSE lineal = {r.rmse_a:.6e} A | RMSE asinh = {r.rmse_log:.6e} | R² lineal = {r.r2:.6f}")
        print(f"  AIC = {r.aic:.3f} | BIC = {r.bic:.3f}")
        for pname, value, err in zip(r.spec.param_names, r.params, r.stderr):
            label, val, sigma = display_param(pname, value, err)
            print(f"  {label:8s} = {val:.6e} ± {sigma:.2e}")

    best = results[0]
    print("\n" + "=" * 78)
    print(f"MEJOR MODELO POR BIC: {best.spec.name}")
    print("Nota: el modelo completo solo debe aceptarse si reduce BIC/AIC y sus parámetros")
    print("son identificables; Rs y Rp no deben sumarse como términos independientes en V(I).")
    print("=" * 78)

    export_theoretical_values(v, results, EXPORT_THEORY_CSV)
    print(f"CSV teórico guardado en: {EXPORT_THEORY_CSV}")

    if GENERAR_GRAFICO:
        make_plot(v, i, results, PLOT_PATH)
        print(f"Gráfico guardado en: {PLOT_PATH}")
    else:
        print("Gráfico omitido porque GENERAR_GRAFICO = False")


if __name__ == "__main__":
    main()
