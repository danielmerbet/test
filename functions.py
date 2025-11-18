"""
GR4J Model with ERA5 Soil Moisture Data Assimilation
Requires: numpy, pandas, scipy, matplotlib
"""

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
from typing import Optional, Dict, Any, Tuple


# -------------------------
# Performance metrics
# -------------------------
def nse(sim: np.ndarray, obs: np.ndarray) -> float:
    """Nash-Sutcliffe efficiency"""
    sim = np.array(sim)
    obs = np.array(obs)
    mask = ~np.isnan(sim) & ~np.isnan(obs)
    if mask.sum() == 0:
        return np.nan
    num = np.sum((obs[mask] - sim[mask]) ** 2)
    den = np.sum((obs[mask] - np.mean(obs[mask])) ** 2)
    return 1 - num / den if den != 0 else np.nan


def rmse(sim: np.ndarray, obs: np.ndarray) -> float:
    mask = ~np.isnan(sim) & ~np.isnan(obs)
    return np.sqrt(np.mean((sim[mask] - obs[mask]) ** 2)) if mask.sum() > 0 else np.nan


def pbias(sim: np.ndarray, obs: np.ndarray) -> float:
    mask = ~np.isnan(sim) & ~np.isnan(obs)
    if mask.sum() == 0 or np.sum(obs[mask]) == 0:
        return np.nan
    return 100.0 * (np.sum(sim[mask]) - np.sum(obs[mask])) / np.sum(obs[mask])


def kge(sim: np.ndarray, obs: np.ndarray) -> float:
    """Kling-Gupta Efficiency (KGE) implementation (version with correlation, alpha, beta)."""
    mask = ~np.isnan(sim) & ~np.isnan(obs)
    if mask.sum() == 0:
        return np.nan
    s = sim[mask]
    o = obs[mask]
    r = np.corrcoef(s, o)[0, 1] if s.size > 1 else 0.0
    alpha = np.std(s) / np.std(o) if np.std(o) != 0 else np.nan
    beta = np.mean(s) / np.mean(o) if np.mean(o) != 0 else np.nan
    # KGE definition (older version): 1 - sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2)
    return 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)


# -------------------------
# Unit hydrograph computation
# -------------------------
def compute_UH(X4: float, type_: int = 1) -> np.ndarray:
    X4 = max(0.5, min(X4, 20.0))
    time_base = X4 if type_ == 1 else 2 * X4
    n = max(1, int(np.ceil(time_base)))
    UH = np.zeros(n)
    for t in range(1, n + 1):
        if t <= time_base:
            val = (t / time_base) ** (5.0 / 2.0)
            if t > 1:
                val -= ((t - 1) / time_base) ** (5.0 / 2.0)
            UH[t - 1] = val
    sum_UH = UH.sum()
    if sum_UH > 0:
        UH = UH / sum_UH
    else:
        UH[0] = 1.0
    return UH


# -------------------------
# GR4J with soil moisture assimilation
# -------------------------
def GR4J_with_SM(
    params: np.ndarray,
    precip: np.ndarray,
    temp: np.ndarray,
    pet: np.ndarray,
    sm_obs: Optional[np.ndarray] = None,
    sm_depths: Optional[np.ndarray] = None,
    warmup: int = 365,
    assimilate: bool = False,
    nudging_weight: float = 1.0,
) -> Dict[str, np.ndarray]:
    """
    params: length-7 array:
      params[0] = X1 (production store max mm)
      params[1] = X2 (GW exchange coefficient mm/day)
      params[2] = X3 (routing store max mm)
      params[3] = X4 (time base UH days)
      params[4] = CTG (snow threshold degC)
      params[5] = Kf (melt factor mm/degC/day)
      params[6] = SM_scale (mm per m3/m3) scaling factor
    """
    precip = np.asarray(precip).astype(float)
    temp = np.asarray(temp).astype(float)
    pet = np.asarray(pet).astype(float)
    n = len(precip)

    # handle missing
    precip[np.isnan(precip)] = 0.0
    if np.isnan(temp).any():
        temp[np.isnan(temp)] = np.nanmean(temp)
    if np.isnan(pet).any():
        pet[np.isnan(pet)] = np.nanmean(pet)

    # allocate
    S = np.zeros(n)  # production store
    R = np.zeros(n)  # routing store
    G = np.zeros(n)  # snow pack mm water equivalent
    eTG = np.zeros(n)
    SM_model = np.zeros(n)

    Pn = np.zeros(n)
    En = np.zeros(n)
    Ps = np.zeros(n)
    Es = np.zeros(n)
    Perc = np.zeros(n)
    Pr = np.zeros(n)
    Q9 = np.zeros(n)
    Q1 = np.zeros(n)
    F = np.zeros(n)
    Qr = np.zeros(n)
    Qd = np.zeros(n)
    Q = np.zeros(n)

    snow_frac = np.zeros(n)
    rain = np.zeros(n)
    snow = np.zeros(n)
    melt = np.zeros(n)

    X4_safe = max(0.5, min(params[3], 20.0))
    UH1 = compute_UH(X4_safe, type_=1)
    UH2 = compute_UH(X4_safe, type_=2)
    len_UH1 = len(UH1)
    len_UH2 = len(UH2)
    store_UH1 = np.zeros(len_UH1)
    store_UH2 = np.zeros(len_UH2)

    # initial conditions
    S[0] = params[0] * 0.5
    R[0] = params[2] * 0.5
    G[0] = 0.0
    eTG[0] = 0.0

    if sm_obs is not None and assimilate:
        # convert observation to store-level (mm)
        S[0] = sm_obs[0] * params[6]

    for t in range(n):
        # snow partition
        if temp[t] <= params[4] - 1.0:
            snow_frac[t] = 1.0
        elif temp[t] >= params[4] + 1.0:
            snow_frac[t] = 0.0
        else:
            snow_frac[t] = (params[4] + 1.0 - temp[t]) / 2.0

        snow[t] = precip[t] * snow_frac[t]
        rain[t] = precip[t] * (1.0 - snow_frac[t])

        # update snow pack
        if t > 0:
            G[t] = G[t - 1] + snow[t]
            eTG[t] = eTG[t - 1]
        else:
            G[t] = snow[t]
            eTG[t] = 0.0

        # snowmelt
        if G[t] > 0.0:
            eTG[t] = eTG[t] + temp[t]
            if eTG[t] > 0.0:
                potential_melt = params[5] * eTG[t]
                melt[t] = min(G[t], potential_melt)
                G[t] = G[t] - melt[t]
                eTG[t] = 0.0
            else:
                melt[t] = 0.0
        else:
            melt[t] = 0.0
            eTG[t] = 0.0

        Pn[t] = rain[t] + melt[t]

        # production store - start from previous
        if t > 0:
            S[t] = S[t - 1]

        # assimilation nudging (note: nudging_weight=1 means full obs)
        if sm_obs is not None and assimilate and not np.isnan(sm_obs[t]):
            S_obs = sm_obs[t] * params[6]
            S[t] = (1.0 - nudging_weight) * S[t] + nudging_weight * S_obs
            S[t] = np.clip(S[t], 0.0, params[0])

        # net evapotranspiration logic
        if Pn[t] >= pet[t]:
            En[t] = 0.0
            scaled_S = S[t] / params[0] if params[0] != 0 else 0.0
            # avoid domain issues with tanh and division by zero:
            denom = 1 + scaled_S * np.tanh(Pn[t] / params[0]) if params[0] != 0 else 1.0
            Ps[t] = params[0] * (1 - scaled_S ** 2) * np.tanh(Pn[t] / params[0]) / denom if params[0] != 0 else 0.0
            Es[t] = 0.0
        else:
            En[t] = pet[t] - Pn[t]
            Ps[t] = 0.0
            scaled_S = S[t] / params[0] if params[0] != 0 else 0.0
            denom = 1 + (1 - scaled_S) * np.tanh(En[t] / params[0]) if params[0] != 0 else 1.0
            Es[t] = S[t] * (2 - scaled_S) * np.tanh(En[t] / params[0]) / denom if params[0] != 0 else 0.0

        # update production store
        S[t] = S[t] - Es[t] + Ps[t]
        S[t] = np.clip(S[t], 0.0, params[0])

        # model soil moisture fraction
        SM_model[t] = S[t] / params[0] if params[0] != 0 else 0.0

        # percolation
        scaled_S = S[t] / params[0] if params[0] != 0 else 0.0
        Perc[t] = S[t] * (1.0 - (1.0 + (scaled_S / 2.25) ** 4) ** (-0.25))
        S[t] = S[t] - Perc[t]
        S[t] = max(0.0, S[t])

        # routing input
        Pr[t] = Perc[t] + (Pn[t] - Ps[t])

        # split for two UH components
        Pr_9 = 0.9 * Pr[t]
        Pr_1 = 0.1 * Pr[t]

        # UH convolution using delay store arrays (shift)
        if len_UH1 > 1:
            store_UH1[1:] = store_UH1[:-1]
        store_UH1[0] = Pr_9

        if len_UH2 > 1:
            store_UH2[1:] = store_UH2[:-1]
        store_UH2[0] = Pr_1

        Q9[t] = np.sum(store_UH1 * UH1)
        Q1[t] = np.sum(store_UH2 * UH2)

        # groundwater exchange & routing store update
        if t > 0:
            R[t] = R[t - 1]
        F[t] = params[1] * (R[t] / params[2]) ** (7.0 / 2.0) if params[2] != 0 else 0.0

        R[t] = max(0.0, R[t] + Q9[t] + F[t])

        Qr[t] = R[t] * (1.0 - (1.0 + (R[t] / params[2]) ** 4) ** (-0.25)) if params[2] != 0 else 0.0
        R[t] = R[t] - Qr[t]
        R[t] = max(0.0, R[t])

        Qd[t] = max(0.0, Q1[t] + F[t])

        Q[t] = Qr[t] + Qd[t]
        Q[t] = max(0.0, Q[t])

    # remove warmup by setting to NaN
    if warmup > 0:
        Q[:warmup] = np.nan

    return {
        "streamflow": Q,
        "production_store": S,
        "routing_store": R,
        "snow_pack": G,
        "snow_thermal_state": eTG,
        "melt": melt,
        "actual_et": Es,
        "percolation": Perc,
        "Qr": Qr,
        "Qd": Qd,
        "rain": rain,
        "snow": snow,
        "Pn": Pn,
        "SM_model": SM_model,
    }


# -------------------------
# Prepare ERA5 SM (weighted root zone)
# -------------------------
def prepare_ERA5_SM(
    sm_0_7: np.ndarray,
    sm_7_28: np.ndarray,
    sm_28_100: np.ndarray,
    sm_100_255: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Weighted average to ~0-50 cm: 7cm + 21cm + 22cm partial of third layer
    Returns volumetric soil moisture (m3/m3)
    """
    depth_1 = 7.0
    depth_2 = 21.0
    depth_3 = 22.0
    total_depth = depth_1 + depth_2 + depth_3

    sm_0_7 = np.asarray(sm_0_7, dtype=float)
    sm_7_28 = np.asarray(sm_7_28, dtype=float)
    sm_28_100 = np.asarray(sm_28_100, dtype=float)

    sm_weighted = (sm_0_7 * depth_1 + sm_7_28 * depth_2 + sm_28_100 * depth_3) / total_depth
    # handle NaNs by replacing with mean
    if np.isnan(sm_weighted).any():
        mean_val = np.nanmean(sm_weighted)
        sm_weighted[np.isnan(sm_weighted)] = mean_val
    print("ERA5 Soil Moisture prepared:")
    print(f"  Mean: {np.nanmean(sm_weighted):.3f} m3/m3")
    print(f"  Range: {np.nanmin(sm_weighted):.3f} - {np.nanmax(sm_weighted):.3f} m3/m3")
    return sm_weighted


# -------------------------
# Objective function for calibration
# -------------------------
def objective_GR4J_SM(
    params_vec: np.ndarray,
    precip: np.ndarray,
    temp: np.ndarray,
    pet: np.ndarray,
    obs_flow: np.ndarray,
    sm_obs: Optional[np.ndarray] = None,
    warmup: int = 365,
    metric: str = "combined",
    assimilate: bool = False,
) -> float:
    # parameter constraints (mimic R checks)
    if params_vec[0] <= 0 or params_vec[2] <= 0 or params_vec[3] <= 0 or params_vec[5] <= 0 or params_vec[6] <= 0:
        return 9999.0

    sim = GR4J_with_SM(params_vec, precip, temp, pet, sm_obs=sm_obs, warmup=warmup, assimilate=assimilate)

    sim_flow = sim["streamflow"]
    # valid indices where both sim and obs exist and obs_flow >= 0
    valid_idx = ~np.isnan(sim_flow) & ~np.isnan(obs_flow) & (obs_flow >= 0)
    sim_valid = sim_flow[valid_idx]
    obs_valid = obs_flow[valid_idx]

    if len(sim_valid) < 30:
        return 9999.0

    if np.any(np.isnan(sim_valid)) or np.any(sim_valid < 0):
        return 9999.0

    if metric == "NSE":
        obj_flow = nse(sim_valid, obs_valid)
    elif metric == "KGE":
        obj_flow = kge(sim_valid, obs_valid)
    else:
        obj_flow = kge(sim_valid, obs_valid)

    # if combined metric with soil moisture
    if (sm_obs is not None) and (metric == "combined"):
        sm_sim = sim["SM_model"]
        valid_sm = (~np.isnan(sm_obs)) & (~np.isnan(sm_sim))
        if np.sum(valid_sm) > 30:
            sm_obs_valid = sm_obs[valid_sm]
            sm_sim_valid = sm_sim[valid_sm]
            sm_corr = np.corrcoef(sm_obs_valid, sm_sim_valid)[0, 1] if sm_obs_valid.size > 1 else 0.0
            sm_rmse = np.sqrt(np.nanmean((sm_obs_valid - sm_sim_valid) ** 2))
            sm_metric = sm_corr - sm_rmse
            obj = 0.8 * obj_flow + 0.2 * sm_metric
        else:
            obj = obj_flow
    else:
        obj = obj_flow

    # negative because we minimize objective in DE; R returned -obj at the end
    return -obj


# -------------------------
# Calibration wrapper using differential evolution
# -------------------------
def calibrate_GR4J_SM(
    precip: np.ndarray,
    temp: np.ndarray,
    pet: np.ndarray,
    obs_flow: np.ndarray,
    sm_obs: Optional[np.ndarray] = None,
    warmup: int = 365,
    max_iter: int = 100,
    metric: str = "KGE",
    assimilate: bool = False,
    lower: Optional[np.ndarray] = None,
    upper: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    # parameter bounds default (7 params)
    if lower is None:
        lower = np.array([10.0, -10.0, 10.0, 0.5, -2.0, 2.0, 100.0])
    if upper is None:
        upper = np.array([1000.0, 10.0, 500.0, 10.0, 2.0, 8.0, 5000.0])

    bounds = [(float(l), float(u)) for l, u in zip(lower, upper)]

    print("=" * 65)
    print("GR4J + SNOW + SOIL MOISTURE CONSTRAINT (Python)")
    print("=" * 65)
    print(f"  X1 (prod store):      {lower[0]:.0f} to {upper[0]:.0f} mm")
    print(f"  X2 (GW exchange):     {lower[1]:.1f} to {upper[1]:.1f} mm/day")
    print(f"  X3 (routing store):   {lower[2]:.0f} to {upper[2]:.0f} mm")
    print(f"  X4 (UH time base):    {lower[3]:.2f} to {upper[3]:.2f} days")
    print(f"  CTG (snow threshold): {lower[4]:.1f} to {upper[4]:.1f} °C")
    print(f"  Kf (melt factor):     {lower[5]:.1f} to {upper[5]:.1f} mm/°C/day")
    print(f"  SM_scale:             {lower[6]:.0f} to {upper[6]:.0f} mm/(m3/m3)\n")

    if sm_obs is not None:
        print("Soil moisture data: AVAILABLE")
        if assimilate:
            print("Data assimilation: ENABLED (nudging factor = 1.0 by default)")
        else:
            print("Data assimilation: DISABLED (used for calibration only)")
    else:
        print("Soil moisture data: NOT PROVIDED")

    print(f"\nCalibration metric: {metric}")
    print(f"Iterations (maxiter): {max_iter}")
    print("Starting calibration ...\n")

    # wrapper for DE
    def obj_wrap(x):
        return objective_GR4J_SM(x, precip, temp, pet, obs_flow, sm_obs=sm_obs, warmup=warmup, metric=("combined" if metric == "combined" else metric), assimilate=assimilate)

    result = differential_evolution(obj_wrap, bounds, maxiter=max_iter, popsize=15, polish=True, disp=True)

    print("\n" + "=" * 65)
    print("Calibration completed!")
    best_params = result.x
    best_val = -result.fun  # because result.fun is negative of objective as we returned -obj
    print(f"Best objective value: {best_val:.4f}\n")
    print("Calibrated parameters:")
    print(f"  X1 (production store):   {best_params[0]:.1f} mm")
    print(f"  X2 (GW exchange):        {best_params[1]:.2f} mm/day")
    print(f"  X3 (routing store):      {best_params[2]:.1f} mm")
    print(f"  X4 (UH time base):       {best_params[3]:.2f} days")
    print(f"  CTG (snow threshold):    {best_params[4]:.2f} °C")
    print(f"  Kf (melt factor):        {best_params[5]:.2f} mm/°C/day")
    print(f"  SM_scale:                {best_params[6]:.1f} mm/(m3/m3)")
    print("=" * 65)
    return {"result": result, "best_params": best_params, "best_val": best_val}


# -------------------------
# Plotting results
# -------------------------
def plot_GR4J_SM_results(
    dates: pd.DatetimeIndex,
    obs_flow: np.ndarray,
    sim_results: Dict[str, np.ndarray],
    sm_obs: Optional[np.ndarray] = None,
    params: Optional[np.ndarray] = None,
    plot_period: Optional[Tuple[pd.Timestamp, pd.Timestamp]] = None,
):
    if plot_period is not None:
        start, end = plot_period
        mask = (dates >= start) & (dates <= end)
        dates = dates[mask]
        obs_flow = obs_flow[mask]
        if sm_obs is not None:
            sm_obs = sm_obs[mask]
        sim_results = {k: v[mask] if len(v) == len(mask) else v for k, v in sim_results.items()}

    n_plots = 5 if sm_obs is not None else 4
    fig, axes = plt.subplots(n_plots, 1, figsize=(16, 3 * n_plots), sharex=True)
    if n_plots == 1:
        axes = [axes]

    # Plot 1: Streamflow + Precip
    ax = axes[0]
    max_flow = np.nanmax(np.concatenate(([np.nanmax(obs_flow)], [np.nanmax(sim_results["streamflow"])])))
    total_precip = sim_results["rain"] + sim_results["snow"]
    max_precip = np.nanmax(total_precip)
    precip_space = max(50.0, max_flow * 0.2)
    y_max = max_flow + precip_space

    ax.plot(dates, obs_flow, label="Observed", color="blue", linewidth=1.5)
    ax.plot(dates, sim_results["streamflow"], label="Simulated", color="red", linewidth=1.0)
    ax.set_ylabel("Streamflow (mm/day)")
    ax.set_ylim(0, y_max)
    ax.set_title("GR4J with Soil Moisture: Streamflow")

    # precipitation as bars mapped to right-side
    if max_precip > 0:
        precip_scale = precip_space / max_precip
        # top of plot is y_max, draw stacks downwards
        top = y_max
        # snow (top) and rain (below)
        ax.bar(dates, sim_results["snow"] * precip_scale, bottom=top - sim_results["snow"] * precip_scale, width=1.0, alpha=0.6, label="Snow")
        ax.bar(dates, sim_results["rain"] * precip_scale, bottom=top - (sim_results["snow"] * precip_scale) - (sim_results["rain"] * precip_scale), width=1.0, alpha=0.7, label="Rain")

        # create secondary axis for precip ticks
        ax2 = ax.twinx()
        ticks = np.linspace(0, max_precip, 5)
        tickpos = y_max - ticks * precip_scale
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(tickpos)
        ax2.set_yticklabels([f"{t:.1f}" for t in ticks])
        ax2.set_ylabel("Precip (mm/day)")

    ax.legend(loc="upper right")

    idx = 1
    # Plot 2: Soil moisture comparison (if available)
    if sm_obs is not None:
        ax_sm = axes[idx]
        ax_sm.plot(dates, sm_obs, label="ERA5 observed", color="darkgreen", linewidth=2)
        ax_sm.plot(dates, sim_results["SM_model"], label="Model", color="brown", linewidth=1.5)
        ax_sm.set_ylabel("Soil moisture (frac)")
        ax_sm.set_title("Soil Moisture: ERA5 vs Model")
        valid = (~np.isnan(sm_obs)) & (~np.isnan(sim_results["SM_model"]))
        if valid.sum() > 10:
            corr = np.corrcoef(sm_obs[valid], sim_results["SM_model"][valid])[0, 1]
            ax_sm.legend(loc="upper right", title=f"r = {corr:.3f}")
        else:
            ax_sm.legend(loc="upper right")
        idx += 1

    # Plot 3: Flow components
    ax3 = axes[idx]
    ax3.plot(dates, sim_results["streamflow"], label="Total", color="black", linewidth=2)
    ax3.plot(dates, sim_results["Qr"], label="Routing", color="blue", linewidth=1)
    ax3.plot(dates, sim_results["Qd"], label="Direct", color="orange", linewidth=1)
    ax3.set_ylabel("Flow (mm/day)")
    ax3.set_title("Flow Components")
    ax3.legend(loc="upper right")
    idx += 1

    # Plot 4: Snow pack
    ax4 = axes[idx]
    ax4.plot(dates, sim_results["snow_pack"], label="Snow pack (mm)", color="darkblue", linewidth=2)
    ax4.set_ylabel("Snow (mm)")
    ax4.set_title("Snow Water Equivalent")
    ax4.legend(loc="upper right")
    idx += 1

    # Plot 5: Production store
    if idx < len(axes):
        ax5 = axes[idx]
        ax5.plot(dates, sim_results["production_store"], label="S", color="brown", linewidth=2)
        if params is not None:
            ax5.axhline(y=params[0], color="red", linestyle="--", label="X1 max")
        ax5.set_ylabel("Store (mm)")
        ax5.set_title("Production Store (S)")
        ax5.legend(loc="upper right")

    plt.xlabel("Date")
    plt.tight_layout()
    plt.show()

# -------------------------
# Calculate metrics and print summary
# -------------------------
def calculate_metrics(obs: np.ndarray, sim: np.ndarray) -> Dict[str, float]:
    mask = ~np.isnan(obs) & ~np.isnan(sim)
    obs_v = obs[mask]
    sim_v = sim[mask]
    metrics = {
        "NSE": nse(sim_v, obs_v),
        "KGE": kge(sim_v, obs_v),
        "RMSE": rmse(sim_v, obs_v),
        "PBIAS": pbias(sim_v, obs_v),
        "R2": np.corrcoef(sim_v, obs_v)[0, 1] ** 2 if obs_v.size > 1 else np.nan,
    }
    print("\n" + "=" * 65)
    print("STREAMFLOW PERFORMANCE:")
    print("=" * 65)
    print(f"NSE:        {metrics['NSE']:.3f}")
    print(f"KGE:        {metrics['KGE']:.3f}")
    print(f"R²:         {metrics['R2']:.3f}")
    print(f"PBIAS:      {metrics['PBIAS']:.1f}%")
    print(f"RMSE:       {metrics['RMSE']:.2f} mm/day")
    print("=" * 65 + "\n")
    return metrics
