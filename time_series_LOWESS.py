import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from statsmodels.nonparametric.smoothers_lowess import lowess

from spn.simulator import Simulator
from spn.models.cxr_txa import build_cxr_txa_net, initial_marking
from spn.visualizer import PetriNetVisualizer


# ============================================================
# Configuration
# ============================================================

SEED = 8
T_MAX = 300.0
MAX_STEPS = 500_000

# Regular sampling interval for the Gillespie trajectory
TIME_STEP = 0.1

# Rolling statistics
ROLLING_WINDOW_HOURS = 100.0

# Ignore the transient before final feature computation
BURN_IN_TIME = 100.0

# LOWESS smoothing parameter:
# larger value -> smoother curve
# smaller value -> follows local fluctuations more closely
LOWESS_FRAC = 0.08

EPSILON = 1e-12


# ============================================================
# Helper functions
# ============================================================

def get_place_value(marking, place_name, place_names):
    """
    Extract a place value from either a dictionary-like
    or vector-like marking.
    """
    try:
        return marking[place_name]
    except (TypeError, KeyError, IndexError):
        index = place_names.index(place_name)
        return marking[index]


def extract_species_trajectory(result, net, species_name):
    """
    Extract one species trajectory from the simulator result.
    """
    place_names = net.place_names()

    return np.asarray(
        [
            get_place_value(marking, species_name, place_names)
            for marking in result.markings
        ],
        dtype=float,
    )


def remove_duplicate_times(times, values):
    """
    If several markings were recorded at the same time,
    retain the last one.
    """
    dataframe = pd.DataFrame(
        {
            "time": times,
            "value": values,
        }
    )

    dataframe = (
        dataframe
        .groupby("time", as_index=False)
        .last()
    )

    return (
        dataframe["time"].to_numpy(dtype=float),
        dataframe["value"].to_numpy(dtype=float),
    )


def resample_stepwise(times, values, regular_times):
    """
    Resample the irregular Gillespie trajectory on a regular grid.

    Gillespie trajectories are piecewise constant, so the value
    at a given time is the most recent marking before that time.
    """
    times, values = remove_duplicate_times(times, values)

    indices = np.searchsorted(
        times,
        regular_times,
        side="right",
    ) - 1

    indices = np.clip(indices, 0, len(values) - 1)

    return values[indices]


def compute_lowess_mean(times, values, frac):
    """
    Estimate the smooth mean/trend curve using LOWESS.
    """
    return lowess(
        endog=values,
        exog=times,
        frac=frac,
        return_sorted=False,
    )


def compute_time_series_features(
    times,
    x1,
    x2,
    rolling_window_hours,
    time_step,
    burn_in_time,
    lowess_frac,
):
    """
    Compute LOWESS means, rolling statistics and final
    stationary features for one simulation run.
    """

    window_size = max(
        2,
        int(round(rolling_window_hours / time_step)),
    )

    dataframe = pd.DataFrame(
        {
            "time": times,
            "x1": x1,
            "x2": x2,
        }
    )

    # --------------------------------------------------------
    # LOWESS smooth mean/trend curves
    # --------------------------------------------------------

    dataframe["x1_lowess_mean"] = compute_lowess_mean(
        times,
        x1,
        lowess_frac,
    )

    dataframe["x2_lowess_mean"] = compute_lowess_mean(
        times,
        x2,
        lowess_frac,
    )

    # --------------------------------------------------------
    # Rolling mean
    # --------------------------------------------------------

    dataframe["x1_rolling_mean"] = (
        dataframe["x1"]
        .rolling(
            window=window_size,
            min_periods=1,
            center=True,
        )
        .mean()
    )

    dataframe["x2_rolling_mean"] = (
        dataframe["x2"]
        .rolling(
            window=window_size,
            min_periods=1,
            center=True,
        )
        .mean()
    )

    # --------------------------------------------------------
    # Rolling standard deviation
    # --------------------------------------------------------

    dataframe["x1_std"] = (
        dataframe["x1"]
        .rolling(
            window=window_size,
            min_periods=2,
            center=True,
        )
        .std(ddof=0)
        .fillna(0.0)
    )

    dataframe["x2_std"] = (
        dataframe["x2"]
        .rolling(
            window=window_size,
            min_periods=2,
            center=True,
        )
        .std(ddof=0)
        .fillna(0.0)
    )

    # --------------------------------------------------------
    # Rolling coefficient of variation
    # --------------------------------------------------------

    dataframe["x1_cv"] = (
        dataframe["x1_std"]
        / (dataframe["x1_rolling_mean"] + EPSILON)
    )

    dataframe["x2_cv"] = (
        dataframe["x2_std"]
        / (dataframe["x2_rolling_mean"] + EPSILON)
    )

    # --------------------------------------------------------
    # Rolling separation score
    # --------------------------------------------------------

    dataframe["separation_score"] = (
        np.abs(
            dataframe["x1_rolling_mean"]
            - dataframe["x2_rolling_mean"]
        )
        / np.sqrt(
            dataframe["x1_std"] ** 2
            + dataframe["x2_std"] ** 2
            + EPSILON
        )
    )

    # --------------------------------------------------------
    # Stationary statistics after burn-in
    # --------------------------------------------------------

    stationary_data = dataframe[
        dataframe["time"] >= burn_in_time
    ]

    if stationary_data.empty:
        raise ValueError(
            f"No observations after burn-in time {burn_in_time}."
        )

    mean_x1 = stationary_data["x1"].mean()
    mean_x2 = stationary_data["x2"].mean()

    std_x1 = stationary_data["x1"].std(ddof=0)
    std_x2 = stationary_data["x2"].std(ddof=0)

    cv_x1 = std_x1 / (mean_x1 + EPSILON)
    cv_x2 = std_x2 / (mean_x2 + EPSILON)

    separation_score = (
        abs(mean_x1 - mean_x2)
        / np.sqrt(
            std_x1 ** 2
            + std_x2 ** 2
            + EPSILON
        )
    )

    # Identify activated and inactive species
    if mean_x1 >= mean_x2:
        activated_species = "x1"
        inactive_species = "x2"

        mean_high = mean_x1
        mean_low = mean_x2

        std_high = std_x1
        std_low = std_x2

        cv_high = cv_x1
        cv_low = cv_x2
    else:
        activated_species = "x2"
        inactive_species = "x1"

        mean_high = mean_x2
        mean_low = mean_x1

        std_high = std_x2
        std_low = std_x1

        cv_high = cv_x2
        cv_low = cv_x1

    summary = {
        "mean_x1": mean_x1,
        "mean_x2": mean_x2,
        "std_x1": std_x1,
        "std_x2": std_x2,
        "cv_x1": cv_x1,
        "cv_x2": cv_x2,
        "separation_score": separation_score,
        "activated_species": activated_species,
        "inactive_species": inactive_species,
        "mean_high": mean_high,
        "mean_low": mean_low,
        "std_high": std_high,
        "std_low": std_low,
        "cv_high": cv_high,
        "cv_low": cv_low,
    }

    return dataframe, summary


# ============================================================
# Plotting functions
# ============================================================

def plot_raw_and_lowess_mean(dataframe):
    """
    Plot the noisy Gillespie trajectories and the LOWESS
    estimated mean/trend curves.
    """
    plt.figure(figsize=(14, 7))

    # Raw trajectories
    plt.plot(
        dataframe["time"],
        dataframe["x1"],
        alpha=0.35,
        linewidth=1.0,
        label="x1 raw",
    )

    plt.plot(
        dataframe["time"],
        dataframe["x2"],
        alpha=0.35,
        linewidth=1.0,
        label="x2 raw",
    )

    # Smooth LOWESS mean curves
    plt.plot(
        dataframe["time"],
        dataframe["x1_lowess_mean"],
        linewidth=3.0,
        label="x1 LOWESS mean",
    )

    plt.plot(
        dataframe["time"],
        dataframe["x2_lowess_mean"],
        linewidth=3.0,
        label="x2 LOWESS mean",
    )

    plt.axvline(
        BURN_IN_TIME,
        linestyle="--",
        linewidth=1.5,
        label="burn-in limit",
    )

    plt.title("Xist dynamics with LOWESS mean curves")
    plt.xlabel("time (h)")
    plt.ylabel("molecule count")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_raw_rolling_and_lowess(dataframe):
    """
    Compare the raw signal, rolling mean and LOWESS mean.
    """
    plt.figure(figsize=(14, 7))

    plt.plot(
        dataframe["time"],
        dataframe["x2"],
        alpha=0.3,
        linewidth=1.0,
        label="x2 raw",
    )

    plt.plot(
        dataframe["time"],
        dataframe["x2_rolling_mean"],
        linewidth=2.0,
        label="x2 rolling mean",
    )

    plt.plot(
        dataframe["time"],
        dataframe["x2_lowess_mean"],
        linewidth=3.0,
        label="x2 LOWESS mean",
    )

    plt.title("Comparison of smoothing methods for x2")
    plt.xlabel("time (h)")
    plt.ylabel("molecule count")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_rolling_coefficient_of_variation(dataframe):
    """
    Plot the rolling coefficient of variation.
    """
    plt.figure(figsize=(14, 6))

    plt.plot(
        dataframe["time"],
        dataframe["x1_cv"],
        label="x1 rolling CV",
    )

    plt.plot(
        dataframe["time"],
        dataframe["x2_cv"],
        label="x2 rolling CV",
    )

    plt.axvline(
        BURN_IN_TIME,
        linestyle="--",
        linewidth=1.5,
        label="burn-in limit",
    )

    plt.title("Rolling coefficient of variation")
    plt.xlabel("time (h)")
    plt.ylabel("coefficient of variation")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_rolling_separation(dataframe):
    """
    Plot the rolling separation score.
    """
    plt.figure(figsize=(14, 6))

    plt.plot(
        dataframe["time"],
        dataframe["separation_score"],
        label="rolling separation score",
    )

    plt.axvline(
        BURN_IN_TIME,
        linestyle="--",
        linewidth=1.5,
        label="burn-in limit",
    )

    plt.title("Separation between x1 and x2")
    plt.xlabel("time (h)")
    plt.ylabel("separation score")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# ============================================================
# Run one Gillespie simulation
# ============================================================

if __name__ == "__main__":

    # Build model
    net = build_cxr_txa_net()
    m0 = initial_marking(net)

    # Run one simulation
    simulator = Simulator(
        net,
        seed=SEED,
    )

    result = simulator.run(
        initial_marking=m0,
        t_max=T_MAX,
        max_steps=MAX_STEPS,
        record_every_step=True,
    )

    # Raw irregular Gillespie times
    event_times = np.asarray(
        result.times,
        dtype=float,
    )

    # Extract x1 and x2
    raw_x1 = extract_species_trajectory(
        result,
        net,
        "x1",
    )

    raw_x2 = extract_species_trajectory(
        result,
        net,
        "x2",
    )

    # Create regular time grid
    regular_times = np.arange(
        0.0,
        T_MAX + TIME_STEP,
        TIME_STEP,
    )

    # Stepwise interpolation of the Gillespie trajectory
    x1 = resample_stepwise(
        event_times,
        raw_x1,
        regular_times,
    )

    x2 = resample_stepwise(
        event_times,
        raw_x2,
        regular_times,
    )

    # Extract time-series features
    feature_dataframe, summary = compute_time_series_features(
        times=regular_times,
        x1=x1,
        x2=x2,
        rolling_window_hours=ROLLING_WINDOW_HOURS,
        time_step=TIME_STEP,
        burn_in_time=BURN_IN_TIME,
        lowess_frac=LOWESS_FRAC,
    )

    # ========================================================
    # Print simulation information
    # ========================================================

    print("Places:")
    print(net.place_names())

    print("\nInitial marking:")
    print(m0)

    print("\nFinal marking:")
    print(result.markings[-1])

    print("\nNumber of fired transitions:")
    print(len(result.fired_transitions))

    print("\nLast 20 fired transitions:")
    print(result.fired_transitions[-20:])

    # ========================================================
    # Print extracted features
    # ========================================================

    print("\n" + "=" * 55)
    print(f"FEATURES AFTER BURN-IN: t >= {BURN_IN_TIME} h")
    print("=" * 55)

    print(f"Mean x1: {summary['mean_x1']:.4f}")
    print(f"Mean x2: {summary['mean_x2']:.4f}")

    print(f"\nStandard deviation x1: {summary['std_x1']:.4f}")
    print(f"Standard deviation x2: {summary['std_x2']:.4f}")

    print(f"\nCoefficient of variation x1: {summary['cv_x1']:.4f}")
    print(f"Coefficient of variation x2: {summary['cv_x2']:.4f}")

    print(
        f"\nSeparation score: "
        f"{summary['separation_score']:.4f}"
    )

    print(
        f"\nActivated species: "
        f"{summary['activated_species']}"
    )

    print(
        f"Inactive species: "
        f"{summary['inactive_species']}"
    )

    print(f"\nMean high state: {summary['mean_high']:.4f}")
    print(f"Mean low state: {summary['mean_low']:.4f}")

    print(f"Standard deviation high state: {summary['std_high']:.4f}")
    print(f"Standard deviation low state: {summary['std_low']:.4f}")

    print(f"CV high state: {summary['cv_high']:.4f}")
    print(f"CV low state: {summary['cv_low']:.4f}")

    # ========================================================
    # Plots
    # ========================================================

    # Raw trajectories plus smooth LOWESS mean curves
    plot_raw_and_lowess_mean(
        feature_dataframe,
    )

    # Compare rolling mean and LOWESS for x2
    plot_raw_rolling_and_lowess(
        feature_dataframe,
    )

    # Rolling CV
    plot_rolling_coefficient_of_variation(
        feature_dataframe,
    )

    # Rolling separation score
    plot_rolling_separation(
        feature_dataframe,
    )

    # Draw final Petri-net marking
    final_marking = result.markings[-1]

    visualizer = PetriNetVisualizer(net)
    visualizer.draw(final_marking)