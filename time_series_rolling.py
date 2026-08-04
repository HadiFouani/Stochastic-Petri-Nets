import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from spn.simulator import Simulator
from spn.models.cxr_txa import build_cxr_txa_net, initial_marking
from spn.visualizer import PetriNetVisualizer


# ============================================================
# Configuration
# ============================================================

SEED = 7
T_MAX = 400.0
MAX_STEPS = 500_000

# Resample the irregular Gillespie trajectory every 0.1 h
TIME_STEP = 0.1

# Rolling statistics over a 10-hour window
ROLLING_WINDOW_HOURS = 299.0

# Ignore the transient before computing final statistics
BURN_IN_TIME = 100.0

EPSILON = 1e-12


# ============================================================
# Helper functions
# ============================================================

def get_place_value(marking, place_name, place_names):
    """
    Extract the value of a place from a marking.

    Supports:
    - dictionary-like markings: marking["x1"]
    - vector/list-like markings: marking[index]
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
    Keep the last recorded marking when multiple events have the same time.
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
    Resample a Gillespie trajectory on a regular time grid.

    Gillespie trajectories are piecewise constant. For each regular time,
    we therefore use the latest marking observed before that time.
    """
    times, values = remove_duplicate_times(times, values)

    indices = np.searchsorted(
        times,
        regular_times,
        side="right",
    ) - 1

    indices = np.clip(indices, 0, len(values) - 1)

    return values[indices]


def compute_time_series_features(
    times,
    x1,
    x2,
    rolling_window_hours,
    time_step,
    burn_in_time,
):
    """
    Compute rolling and stationary features for one simulation run.
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

    # Rolling mean
    dataframe["x1_mean"] = (
        dataframe["x1"]
        .rolling(
            window=window_size,
            min_periods=1,
        )
        .mean()
    )

    dataframe["x2_mean"] = (
        dataframe["x2"]
        .rolling(
            window=window_size,
            min_periods=1,
        )
        .mean()
    )

    # Rolling standard deviation
    dataframe["x1_std"] = (
        dataframe["x1"]
        .rolling(
            window=window_size,
            min_periods=2,
        )
        .std(ddof=0)
        .fillna(0.0)
    )

    dataframe["x2_std"] = (
        dataframe["x2"]
        .rolling(
            window=window_size,
            min_periods=2,
        )
        .std(ddof=0)
        .fillna(0.0)
    )

    # Rolling coefficient of variation
    dataframe["x1_cv"] = (
        dataframe["x1_std"]
        / (dataframe["x1_mean"] + EPSILON)
    )

    dataframe["x2_cv"] = (
        dataframe["x2_std"]
        / (dataframe["x2_mean"] + EPSILON)
    )

    # Rolling separation score
    dataframe["separation_score"] = (
        np.abs(dataframe["x1_mean"] - dataframe["x2_mean"])
        / np.sqrt(
            dataframe["x1_std"] ** 2
            + dataframe["x2_std"] ** 2
            + EPSILON
        )
    )

    # Long-term statistics after burn-in
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

    # Identify the activated and inactive allele
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


def plot_raw_and_rolling_mean(
    dataframe,
    rolling_window_hours,
):
    """
    Plot raw x1/x2 trajectories and their rolling means.
    """
    plt.figure(figsize=(14, 7))

    plt.plot(
        dataframe["time"],
        dataframe["x1"],
        alpha=0.35,
        linewidth=1.0,
        label="x1",
    )

    plt.plot(
        dataframe["time"],
        dataframe["x2"],
        alpha=0.35,
        linewidth=1.0,
        label="x2",
    )

    plt.plot(
        dataframe["time"],
        dataframe["x1_mean"],
        linewidth=2.5,
        label=f"x1 rolling mean ",
    )

    plt.plot(
        dataframe["time"],
        dataframe["x2_mean"],
        linewidth=2.5,
        label=f"x2 rolling mean ",
    )

    plt.axvline(
        BURN_IN_TIME,
        linestyle="--",
        linewidth=1.5,
        label="burn-in limit",
    )

    plt.title("Xist trajectories and rolling means")
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

    # Build model and initial marking
    net = build_cxr_txa_net()
    m0 = initial_marking(net)

    # Run simulation
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

    # Extract raw event times and species trajectories
    event_times = np.asarray(
        result.times,
        dtype=float,
    )

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

    # Gillespie event times are irregular, so construct a regular grid
    regular_times = np.arange(
        0.0,
        T_MAX + TIME_STEP,
        TIME_STEP,
    )

    # Piecewise-constant resampling
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

    # Compute rolling and stationary features
    feature_dataframe, summary = compute_time_series_features(
        times=regular_times,
        x1=x1,
        x2=x2,
        rolling_window_hours=ROLLING_WINDOW_HOURS,
        time_step=TIME_STEP,
        burn_in_time=BURN_IN_TIME,
    )

    # Print simulation information
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

    # Print extracted features
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

    print(f"CV high state: {summary['cv_high']:.4f}")
    print(f"CV low state: {summary['cv_low']:.4f}")

    # Plot the raw trajectories and rolling means
    plot_raw_and_rolling_mean(
        feature_dataframe,
        ROLLING_WINDOW_HOURS,
    )

    # Plot rolling CV
    plot_rolling_coefficient_of_variation(
        feature_dataframe,
    )

    # Plot rolling separation score
    plot_rolling_separation(
        feature_dataframe,
    )

    # Draw final Petri-net marking
    final_marking = result.markings[-1]

    visualizer = PetriNetVisualizer(net)
    visualizer.draw(final_marking)