"""MSC feedback model.

Document-consistent Material-Selective Capture (MSC) simulator with an optional
slow Hill-feedback material tag. This is intended for Google Colab runs.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
import math
import json
import numpy as np
import pandas as pd

try:
    import yaml
except Exception:  # pragma: no cover
    yaml = None


MODEL_VARIANTS: Dict[str, Dict[str, Any]] = {
    "full_msc_feedback": {
        "label": "Full MSC + feedback",
        "material": True,
        "feedback": True,
        "mode": "correct",
        "axes": "primary3",
    },
    "msc_capture_only": {
        "label": "MSC capture only",
        "material": True,
        "feedback": False,
        "mode": "correct",
        "axes": "primary3",
    },
    "projected1d_feedback": {
        "label": "Projected 1D + feedback",
        "material": True,
        "feedback": True,
        "mode": "correct",
        "axes": "axis0",
    },
    "shuffled_material_feedback": {
        "label": "Shuffled material + feedback",
        "material": True,
        "feedback": True,
        "mode": "shuffled",
        "axes": "primary3",
    },
    "random_material_feedback": {
        "label": "Random material + feedback",
        "material": True,
        "feedback": True,
        "mode": "random",
        "axes": "primary3",
    },
    "tag_only_feedback": {
        "label": "Activity/tag + feedback only",
        "material": False,
        "feedback": True,
        "mode": "none",
        "axes": "primary3",
    },
}

DEFAULT_RANGES: Dict[str, Tuple[float, float]] = {
    "L": (80.0, 180.0),
    "Nsyn": (18, 35),
    "granule_density": (0.60, 1.30),
    "v": (0.05, 0.30),
    "D": (0.005, 0.15),
    "p_switch": (0.002, 0.04),
    "epsilon": (0.20, 0.80),
    "h0": (5e-4, 2.5e-3),
    "beta_T": (2.0, 3.8),
    "beta_M": (1.0, 3.0),
    "ell_h": (1.0, 2.2),
    "ell_q_ratio": (0.35, 0.75),
    "alphaS": (0.025, 0.080),
    "q_min": (0.10, 0.25),
    "tau_T": (900.0, 3600.0),
    "tau_S": (1200.0, 7200.0),
    "delta_translation": (300.0, 900.0),
    "M_bas": (0.02, 0.12),
    "M_max": (0.85, 1.00),
    "k_M": (1.5e-4, 1.2e-3),
    "alphaM": (0.0015, 0.010),
    "K_M": (0.40, 0.65),
    "n_hill": (3, 6),
    "strong_fraction": (0.35, 0.60),
    "strong_sigma": (0.18, 0.42),
    "weak_shift": (1.2, 2.4),
    "weak_sigma": (0.65, 1.10),
}

DEFAULT_FIXED: Dict[str, Any] = {
    "D_material": 4,
    "dt": 2.0,
    "T_total": 3600.0,
    "activation_time": 120.0,
    "activation_amp": 1.0,
    "late_window": 600.0,
    "S_threshold": 0.60,
    "M_threshold": 0.50,
    "cargo_ell": 1.0,
    "max_obs_per_run": 5000,
}


def _as_range(value: Any) -> Optional[Tuple[float, float]]:
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return float(value[0]), float(value[1])
    return None


def load_config(path: Optional[str | Path]) -> Dict[str, Any]:
    if path is None:
        return {}
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() in {".yaml", ".yml"}:
        if yaml is None:
            raise RuntimeError("pyyaml is required to read YAML config files")
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
    elif path.suffix.lower() == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
    else:
        raise ValueError(f"Unsupported config suffix: {path.suffix}")
    return data


def sample_params(seed: int, config: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    config = config or {}
    rng = np.random.default_rng(seed)
    ranges = dict(DEFAULT_RANGES)
    fixed = dict(DEFAULT_FIXED)

    for key, value in (config.get("parameter_ranges") or {}).items():
        rr = _as_range(value)
        if rr is not None:
            ranges[key] = rr
    for key, value in (config.get("fixed") or {}).items():
        fixed[key] = value

    p: Dict[str, Any] = dict(fixed)
    for key, (lo, hi) in ranges.items():
        if key in {"Nsyn", "n_hill"}:
            p[key] = int(rng.integers(int(round(lo)), int(round(hi)) + 1))
        else:
            p[key] = float(rng.uniform(float(lo), float(hi)))

    p["Ngran"] = int(max(10, round(p["L"] * p["granule_density"])))
    p["ell_q"] = max(0.05, p["ell_h"] * p["ell_q_ratio"])
    p["seed"] = int(seed)
    return p


def axes_from_variant(variant: Dict[str, Any], D: int) -> np.ndarray:
    axes = variant.get("axes", "primary3")
    if axes == "axis0":
        return np.array([0], dtype=int)
    if axes == "primary3":
        return np.arange(min(3, D), dtype=int)
    if axes == "all":
        return np.arange(D, dtype=int)
    if isinstance(axes, (list, tuple)):
        return np.array(axes, dtype=int)
    return np.arange(min(3, D), dtype=int)


def material_distance2(u: np.ndarray, theta: np.ndarray, axes: np.ndarray) -> float:
    diff = u[axes] - theta[axes]
    return float(np.mean(diff * diff))


def make_material_pool(p: Dict[str, Any], rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    D = int(p["D_material"])
    Ngran = int(p["Ngran"])
    n_strong = int(np.clip(round(Ngran * p["strong_fraction"]), 1, Ngran - 1))
    u_strong = rng.normal(0.0, p["strong_sigma"], size=(n_strong, D))
    signs = rng.choice([-1.0, 1.0], size=(Ngran - n_strong, D))
    u_weak = signs * rng.normal(p["weak_shift"], p["weak_sigma"], size=(Ngran - n_strong, D))
    u = np.vstack([u_strong, u_weak])
    perm = rng.permutation(Ngran)
    u = u[perm]
    d2_to_active = np.mean(u[:, : min(3, D)] ** 2, axis=1)
    q_g = np.exp(-0.5 * d2_to_active / (float(p["cargo_ell"]) ** 2))
    strong = q_g >= np.quantile(q_g, 1.0 - float(p["strong_fraction"]))
    return u, q_g.astype(float), strong.astype(bool)


def reservoir_add(rows: List[Dict[str, Any]], row: Dict[str, Any], seen: int, max_rows: int, rng: np.random.Generator) -> int:
    seen += 1
    if max_rows <= 0:
        return seen
    if len(rows) < max_rows:
        rows.append(row)
    else:
        j = int(rng.integers(0, seen))
        if j < max_rows:
            rows[j] = row
    return seen


def simulate_one(
    params: Dict[str, Any],
    seed: int,
    variant_name: str = "full_msc_feedback",
    active_synapses: Optional[Sequence[int]] = None,
    activation_amp: Optional[float] = None,
    initial_S: Optional[Dict[int, float]] = None,
    initial_M: Optional[Dict[int, float]] = None,
    record_obs: bool = False,
    record_trace: bool = False,
) -> Dict[str, Any]:
    if variant_name not in MODEL_VARIANTS:
        raise ValueError(f"Unknown variant {variant_name}. Options: {list(MODEL_VARIANTS)}")
    variant = MODEL_VARIANTS[variant_name]
    p = dict(params)
    rng = np.random.default_rng(int(seed))

    Nsyn = int(p["Nsyn"])
    Ngran = int(p["Ngran"])
    L = float(p["L"])
    Dm = int(p["D_material"])
    dt = float(p["dt"])
    n_steps = int(round(float(p["T_total"]) / dt))
    activation_step = int(round(float(p["activation_time"]) / dt))
    delay_steps = max(1, int(round(float(p["delta_translation"]) / dt)))
    late_steps = max(1, int(round(float(p["late_window"]) / dt)))
    max_obs = int(p.get("max_obs_per_run", 5000))
    activation_amp = float(p["activation_amp"] if activation_amp is None else activation_amp)

    xsyn = np.linspace(L / (Nsyn + 1), L - L / (Nsyn + 1), Nsyn)
    if active_synapses is None:
        active_synapses = [Nsyn // 2]
    active_synapses = [int(x) for x in active_synapses]
    active_mask = np.zeros(Nsyn, dtype=bool)
    active_mask[active_synapses] = True

    theta = rng.normal(0.0, 1.0, size=(Nsyn, Dm))
    for i in active_synapses:
        theta[i, :] = 0.0

    u_true, q_g, strong = make_material_pool(p, rng)
    if variant["mode"] == "shuffled":
        u_used = u_true[rng.permutation(Ngran)]
    elif variant["mode"] == "random":
        u_used = rng.normal(0.0, 1.25, size=u_true.shape)
    else:
        u_used = u_true.copy()

    gpos = rng.uniform(0.0, L, size=Ngran)
    gdir = rng.choice([-1.0, 1.0], size=Ngran)

    S = np.zeros(Nsyn, dtype=float)
    M = np.full(Nsyn, float(p["M_bas"]), dtype=float)
    Ttag = np.zeros(Nsyn, dtype=float)
    if initial_S:
        for k, v in initial_S.items():
            S[int(k)] = float(v)
    if initial_M:
        for k, v in initial_M.items():
            M[int(k)] = float(v)

    pending = np.zeros((n_steps + delay_steps + 2, Nsyn), dtype=float)
    S_hist = np.zeros((n_steps, Nsyn), dtype=np.float32) if (record_trace or True) else None
    M_hist = np.zeros((n_steps, Nsyn), dtype=np.float32) if (record_trace or True) else None
    T_active_hist = np.zeros(n_steps, dtype=np.float32) if record_trace else None

    axes = axes_from_variant(variant, Dm)
    dock_rows: List[Dict[str, Any]] = []
    obs_rows: List[Dict[str, Any]] = []
    obs_seen = 0

    feedback_on = bool(variant["feedback"])
    material_on = bool(variant["material"])
    Kn = float(p["K_M"]) ** int(p["n_hill"])

    for step in range(n_steps):
        if step == activation_step:
            for i in active_synapses:
                Ttag[i] = min(1.0, Ttag[i] + activation_amp)

        S += pending[step]
        S += (-S / float(p["tau_S"])) * dt
        S = np.clip(S, 0.0, 1.0)

        if feedback_on:
            Sn = S ** int(p["n_hill"])
            H = Sn / (Kn + Sn + 1e-12)
            M += (-float(p["k_M"]) * (M - float(p["M_bas"])) + float(p["alphaM"]) * H * (float(p["M_max"]) - M)) * dt
            M = np.clip(M, 0.0, 1.0)

        Ttag += (-Ttag / float(p["tau_T"])) * dt
        Ttag = np.clip(Ttag, 0.0, 1.0)

        gpos += gdir * float(p["v"]) * dt + np.sqrt(2.0 * float(p["D"]) * dt) * rng.standard_normal(Ngran)
        flips = rng.random(Ngran) < float(p["p_switch"])
        gdir[flips] *= -1.0
        below = gpos < 0.0
        above = gpos > L
        gpos[below] = -gpos[below]
        gpos[above] = 2.0 * L - gpos[above]
        gdir[below | above] *= -1.0

        for g in range(Ngran):
            near = np.flatnonzero(np.abs(xsyn - gpos[g]) <= float(p["epsilon"]))
            if near.size == 0:
                continue
            for i in near:
                d2 = material_distance2(u_used[g], theta[i], axes)
                if material_on:
                    phi_h = math.exp(-0.5 * d2 / (float(p["ell_h"]) ** 2))
                    phi_q = math.exp(-0.5 * d2 / (float(p["ell_q"]) ** 2))
                else:
                    phi_h = 1.0
                    phi_q = 1.0
                exponent = float(p["beta_T"]) * Ttag[i]
                if feedback_on:
                    exponent += float(p["beta_M"]) * M[i]
                hazard = float(p["h0"]) * math.exp(min(20.0, exponent)) * phi_h
                pdock = 1.0 - math.exp(-hazard * dt)
                dock = rng.random() < pdock

                if record_obs:
                    obs_seen = reservoir_add(
                        obs_rows,
                        {
                            "seed": int(seed),
                            "variant": variant_name,
                            "T": float(Ttag[i]),
                            "M": float(M[i]),
                            "d2": float(d2),
                            "active_synapse": int(active_mask[i]),
                            "dock": int(dock),
                        },
                        obs_seen,
                        max_obs,
                        rng,
                    )

                if dock:
                    qig = float(q_g[g] * phi_q)
                    productive = qig >= float(p["q_min"])
                    if productive:
                        deliver = min(step + delay_steps, pending.shape[0] - 1)
                        pending[deliver, i] += float(p["alphaS"]) * qig
                    dock_rows.append(
                        {
                            "seed": int(seed),
                            "variant": variant_name,
                            "t": step * dt,
                            "synapse": int(i),
                            "granule": int(g),
                            "active_synapse": int(active_mask[i]),
                            "d2": float(d2),
                            "q_g": float(q_g[g]),
                            "q_ig": float(qig),
                            "strong": int(strong[g]),
                            "productive": int(productive),
                            "T": float(Ttag[i]),
                            "M": float(M[i]),
                        }
                    )
        S_hist[step] = S
        M_hist[step] = M
        if record_trace:
            T_active_hist[step] = float(Ttag[active_mask].mean())

    late_slice = slice(max(0, n_steps - late_steps), n_steps)
    active_S_late = float(S_hist[late_slice][:, active_mask].mean())
    control_late = S_hist[late_slice][:, ~active_mask]
    control_max_S_late = float(control_late.mean(axis=0).max()) if control_late.size else 0.0
    active_M_late = float(M_hist[late_slice][:, active_mask].mean())
    control_M_late = M_hist[late_slice][:, ~active_mask]
    control_max_M_late = float(control_M_late.mean(axis=0).max()) if control_M_late.size else 0.0

    dock_df = pd.DataFrame(dock_rows)
    pool_strong_fraction = float(np.mean(strong))
    pool_productive_weight_mean = float(np.mean(q_g))
    if len(dock_df):
        active_docks_df = dock_df[dock_df["active_synapse"] == 1]
        active_docks = int(len(active_docks_df))
        n_productive = int(dock_df["productive"].sum())
        active_productive_docks = int(active_docks_df["productive"].sum()) if active_docks else 0
        mean_material_distance_docked = float(np.sqrt(dock_df["d2"].astype(float)).mean())
    else:
        active_docks_df = pd.DataFrame()
        active_docks = 0
        n_productive = 0
        active_productive_docks = 0
        mean_material_distance_docked = np.nan

    if active_docks > 0:
        active_strong_dock_fraction = float(active_docks_df["strong"].mean())
        active_productive_weight_mean = float(active_docks_df["q_ig"].mean())
    else:
        active_strong_dock_fraction = np.nan
        active_productive_weight_mean = np.nan

    metrics = {
        "seed": int(params.get("seed", seed)),
        "sim_seed": int(seed),
        "variant": variant_name,
        "label": variant["label"],
        "active_S_late": active_S_late,
        "control_max_S_late": control_max_S_late,
        "continuous_specificity": active_S_late - control_max_S_late,
        "active_M_late": active_M_late,
        "control_max_M_late": control_max_M_late,
        "material_memory_specificity": active_M_late - control_max_M_late,
        "active_consolidation": float(active_S_late >= float(p["S_threshold"])),
        "active_material_memory": float(active_M_late >= float(p["M_threshold"])),
        "crosstalk_rate": float(control_max_S_late >= float(p["S_threshold"])),
        "material_crosstalk_rate": float(control_max_M_late >= float(p["M_threshold"])),
        "n_docks": int(len(dock_df)),
        "n_productive_docks": n_productive,
        "active_docks": active_docks,
        "active_productive_docks": active_productive_docks,
        "pool_strong_fraction": pool_strong_fraction,
        "active_strong_dock_fraction": active_strong_dock_fraction,
        "active_strong_enrichment": active_strong_dock_fraction / pool_strong_fraction if active_docks > 0 and pool_strong_fraction > 0 else np.nan,
        "pool_productive_weight_mean": pool_productive_weight_mean,
        "active_productive_weight_mean": active_productive_weight_mean,
        "active_productive_enrichment": active_productive_weight_mean / pool_productive_weight_mean if active_docks > 0 and pool_productive_weight_mean > 0 else np.nan,
        "mean_material_distance_docked": mean_material_distance_docked,
    }
    for key in ["Nsyn", "Ngran", "L", "granule_density", "v", "D", "epsilon", "h0", "beta_T", "beta_M", "ell_h", "ell_q", "alphaS", "q_min", "tau_T", "tau_S", "delta_translation", "k_M", "alphaM", "K_M", "n_hill"]:
        metrics[key] = p.get(key, np.nan)

    return {
        "metrics": metrics,
        "dock_df": dock_df,
        "obs_df": pd.DataFrame(obs_rows),
        "S_history": S_hist if record_trace else None,
        "M_history": M_hist if record_trace else None,
        "T_active_history": T_active_hist if record_trace else None,
        "xsyn": xsyn,
        "params": p,
    }


def summarize_metrics(metrics_df: pd.DataFrame) -> pd.DataFrame:
    numeric_cols = [c for c in metrics_df.columns if c not in {"variant", "label"} and pd.api.types.is_numeric_dtype(metrics_df[c])]
    rows = []
    for variant, sub in metrics_df.groupby("variant", sort=False):
        row: Dict[str, Any] = {"variant": variant, "label": sub["label"].iloc[0], "n": int(len(sub))}
        for col in numeric_cols:
            vals = sub[col].replace([np.inf, -np.inf], np.nan).dropna()
            if len(vals) == 0:
                continue
            row[f"{col}_mean"] = float(vals.mean())
            row[f"{col}_sem"] = float(vals.std(ddof=1) / math.sqrt(len(vals))) if len(vals) > 1 else 0.0
            row[f"{col}_median"] = float(vals.median())
        rows.append(row)
    return pd.DataFrame(rows)


def paired_effects(metrics_df: pd.DataFrame, reference: str = "full_msc_feedback") -> pd.DataFrame:
    compare_metrics = [
        "active_strong_enrichment",
        "active_productive_enrichment",
        "continuous_specificity",
        "material_memory_specificity",
        "crosstalk_rate",
        "material_crosstalk_rate",
        "n_docks",
        "active_docks",
    ]
    rows = []
    for metric in compare_metrics:
        pivot = metrics_df.pivot_table(index="seed", columns="variant", values=metric, aggfunc="mean")
        if reference not in pivot.columns:
            continue
        for variant in pivot.columns:
            if variant == reference:
                continue
            diff = (pivot[reference] - pivot[variant]).replace([np.inf, -np.inf], np.nan).dropna()
            if len(diff) == 0:
                continue
            rows.append({
                "metric": metric,
                "reference": reference,
                "comparison": variant,
                "n_pairs": int(len(diff)),
                "mean_difference": float(diff.mean()),
                "sem_difference": float(diff.std(ddof=1) / math.sqrt(len(diff))) if len(diff) > 1 else 0.0,
                "median_difference": float(diff.median()),
                "fraction_reference_greater": float((diff > 0).mean()),
            })
    return pd.DataFrame(rows)


def fit_binary_observation_model(obs_df: pd.DataFrame) -> pd.DataFrame:
    if obs_df.empty or obs_df["dock"].nunique() < 2:
        return pd.DataFrame([{"term": "fit_failed", "estimate": np.nan, "note": "empty or single-class observation data"}])
    try:
        import statsmodels.api as sm
        X = obs_df[["T", "M", "d2", "active_synapse"]].astype(float)
        X = sm.add_constant(X, has_constant="add")
        y = obs_df["dock"].astype(float)
        fit = sm.GLM(y, X, family=sm.families.Binomial()).fit(maxiter=100, disp=False)
        names = ["Intercept" if n == "const" else n for n in fit.params.index]
        expected = {"T": "+", "M": "+", "d2": "-", "active_synapse": "weak/absorbed by T"}
        return pd.DataFrame({
            "term": names,
            "estimate": fit.params.values,
            "std_error": fit.bse.values,
            "z": fit.tvalues.values,
            "pvalue": fit.pvalues.values,
            "expected_sign": [expected.get(n, "") for n in names],
            "n_rows": int(len(obs_df)),
            "n_docks": int(obs_df["dock"].sum()),
        })
    except Exception as e:
        return pd.DataFrame([{"term": "fit_failed", "estimate": np.nan, "note": str(e), "n_rows": len(obs_df), "n_docks": int(obs_df["dock"].sum())}])


def equation_mapping() -> pd.DataFrame:
    rows = [
        ("Eq. 18", "Encounter gate G_ig(t)", "simulate_one(): near = abs(xsyn - gpos[g]) <= epsilon", "implemented"),
        ("Eq. 23", "Docking compatibility kernel Phi_h", "phi_h = exp(-0.5*d2/ell_h^2)", "implemented"),
        ("Eq. 26 / 77", "Continuous-time MSC docking hazard", "hazard = h0 * exp(beta_T*T + beta_M*M) * Phi_h", "implemented; feedback factor optional"),
        ("Eq. 32 / 78", "Hazard-to-probability conversion", "Pdock = 1 - exp(-hazard*dt)", "implemented"),
        ("Eq. 33 / 79", "Productive payload", "q_ig = q_g * exp(-0.5*d2/ell_q^2)", "implemented"),
        ("Eq. 34 / 80", "Delayed reinforcement", "pending[step + delay_steps] += alphaS*q_ig", "implemented"),
        ("Eq. 35", "Fast activity tag decay", "Ttag += (-Ttag/tau_T)*dt", "implemented"),
        ("Eq. 38/39", "Synaptic state dynamics", "S += -S/tau_S*dt + delayed inputs", "implemented"),
        ("Eq. 40/41", "Slow material feedback tag", "M += [-k_M(M-M_bas)+alphaM*H(S)(M_max-M)]dt", "implemented"),
        ("Eq. 64", "Binary observation model", "dock ~ T + M + d2 + active_synapse", "implemented as analysis"),
        ("Section 14.2", "Avoid circularity", "primary3 axes exclude outcome-derived dwell axis", "implemented by default"),
    ]
    return pd.DataFrame(rows, columns=["document_equation", "concept", "code_location", "status"])
