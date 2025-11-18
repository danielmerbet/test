
import os
from matplotlib.backends.backend_pdf import PdfPages

os.chdir("/home/dmercado/Documents/NEREIDA/test/")
from functions import *
if __name__ == "__main__":
    # EDIT THESE PATHS & BASIN NAMES TO YOUR LOCAL FILES BEFORE RUNNING
    dir_path = "/home/dmercado/Documents/NEREIDA/test/"
    basin = "ter"  # or "tordera"

    # Load main data (assumes same structure as your R CSV)
    df = pd.read_csv(f"{dir_path}/{basin}/data/water/iwrm/data_iwrm.csv", parse_dates=["date"])
    area = 1661.609  # basin area used in R script
    data = pd.DataFrame({
        "date": df["date"],
        "precip": df["P_Sau"],
        "temp": df["T_Sau"],
        "evap": df["E_Sau"],
        "obs_flow": df["Qin_calc_sau"] * 1000.0 * 86400.0 / (area * 1e6),
    })

    dates = pd.to_datetime(data["date"])
    precip = data["precip"].values
    temp = data["temp"].values
    pet = data["evap"].values
    obs_flow = data["obs_flow"].values

    # soil file (assumes similar columns)
    soil = pd.read_csv(f"{dir_path}/{basin}/data/meteo/meteo-soil_dailydata_predam.csv", parse_dates=["date"])
    soil_data = pd.merge(data, soil, on="date", how="inner")

    sm_0_7 = soil_data["soil_moisture_0_to_7cm"].values
    sm_7_28 = soil_data["soil_moisture_7_to_28cm"].values
    sm_28_100 = soil_data["soil_moisture_28_to_100cm"].values
    # sm_100_255 available but not used in weighted function
    sm_combined = prepare_ERA5_SM(sm_0_7, sm_7_28, sm_28_100)

    warmup = 365
    calibrate = False
    nudge = 100
    subbasin = 'ter_sau'

    # load previously saved best_params (numpy) OR calibrate

    model_dir = os.path.join(dir_path, basin, "models", "river", "gr4j_sm")
    os.makedirs(model_dir, exist_ok=True)
    params_file = os.path.join(model_dir, f"par_nudg{nudge}_{subbasin}.npy")

    if calibrate:
        res = calibrate_GR4J_SM(
            precip,
            temp,
            pet,
            obs_flow,
            sm_obs=sm_combined,
            warmup=warmup,
            max_iter=50,  # set lower/higher depending on time
            metric="combined",
            assimilate=True,
        )
        best_params = res["best_params"]
        np.save(params_file, best_params)
    else:
        # load parameters saved previously (ensure file exists)
        if os.path.exists(params_file):
            best_params = np.load(params_file)
        #else:
            # fallback: some example reasonable params (replace with your own)
         #   best_params = np.array([300.0, 0.5, 50.0, 3.0, 0.0, 4.0, 200.0])

    sim = GR4J_with_SM(best_params, precip, temp, pet, sm_obs=sm_combined, assimilate=True)

    # Plot and metrics
    plt.figure(figsize=(3, 6))
    plot_GR4J_SM_results(dates, obs_flow, sim, sm_obs=sm_combined, params=best_params)
    plt.savefig(model_dir + '/sau_py.pdf', bbox_inches='tight', dpi=300)
    _metrics = calculate_metrics(obs_flow, sim["streamflow"])

    # Simple volume balance example (Method 1)
    bal_run = pd.read_csv(f"{dir_path}/{basin}/data/water/iwrm/data_iwrm.csv", parse_dates=["date"])
    # drop warmup rows to keep in line with R script
    bal_run = bal_run.iloc[warmup:].reset_index(drop=True)

    # Convert Qsim (mm/day) to m3/s or m3/day depending on usage as in R:
    # In R: Qsim_sau <- sim$streamflow*(area*1e6)/(1000)  (mm -> m3) then uses in volume balance
    Qsim_sau = sim["streamflow"] * (area * 1e6) / 1000.0
    Qsim_sau = Qsim_sau[warmup:]  # match R logic

    # Vini: initial reservoir volume from file (assumes V_sau column is in millions m3 in R)
    Vini = bal_run.loc[0, "V_sau"] * 1e6 if "V_sau" in bal_run.columns else 0.0
    Qout_sau = bal_run["Qout_calc_sau"].values * 86400.0 if "Qout_calc_sau" in bal_run.columns else np.zeros_like(Qsim_sau)
    Vtemp = Vini
    Vall = [Vini]
    capacity = 166e6  # as in R

    for v, qsim in enumerate(Qsim_sau[1:], start=1):
        Qout_val = Qout_sau[v] if v < len(Qout_sau) else 0.0
        Vtemp = max(0.0, Vtemp + qsim - Qout_val)
        Vtemp = min(capacity, Vtemp)
        Vall.append(Vtemp)

    # quick plot compare
    p_s = plt.figure()
    plt.plot(pd.to_datetime(bal_run["date"]), Vall, label="Simulated Volume")
    if "V_sau" in bal_run.columns:
        plt.plot(pd.to_datetime(bal_run["date"]), bal_run["V_sau"] * 1e6, label="Observed V_sau", color="blue")
    plt.legend()
    plt.show()
    p_s.savefig("foo.pdf", bbox_inches='tight')
