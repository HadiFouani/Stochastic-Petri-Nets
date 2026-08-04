import matplotlib.pyplot as plt
import numpy as np

from spn.simulator import Simulator
from spn.models.cxr_txa import build_cxr_txa_net, initial_marking, p21


NUM_RUNS = 10
T_MAX = 300.0
MAX_STEPS = 50000

# Analysis settings
XIST_HIGH_FRACTION = 0.2
XIST_HIGH_THRESHOLD = XIST_HIGH_FRACTION * p21
INITIAL_SLOPE_WINDOW = 30.0  # hours


def time_weighted_mean(times, values, start=0.0, end=T_MAX):
    """Mean of a piecewise-constant Gillespie trajectory on [start, end]."""
    boundaries = np.concatenate(([start], times[(times > start) & (times < end)], [end]))
    indices = np.searchsorted(times, boundaries[:-1], side="right") - 1
    indices = np.clip(indices, 0, len(values) - 1)
    durations = np.diff(boundaries)
    return float(np.sum(values[indices] * durations) / (end - start))


def initial_slope(times, values, window=INITIAL_SLOPE_WINDOW):
    """Linear slope fitted to regularly sampled values at the start of a run."""
    end = min(float(times[-1]), window)
    if end <= 0.0:
        return 0.0

    sample_times = np.linspace(0.0, end, 101)
    indices = np.searchsorted(times, sample_times, side="right") - 1
    indices = np.clip(indices, 0, len(values) - 1)
    sampled_values = values[indices]
    return float(np.polyfit(sample_times, sampled_values, 1)[0])


def classify_choice(x1_mean, x2_mean):
    x1_high = x1_mean >= XIST_HIGH_THRESHOLD
    x2_high = x2_mean >= XIST_HIGH_THRESHOLD

    if x1_high and not x2_high:
        return "x1"
    if x2_high and not x1_high:
        return "x2"
    if x1_high and x2_high:
        return "both"
    return "neither"


if __name__ == "__main__":
    net = build_cxr_txa_net()
    m0 = initial_marking(net)
    names = net.place_names()
    x1_index = names.index("x1")
    x2_index = names.index("x2")

    summaries = []
    for run in range(NUM_RUNS):
        simulator = Simulator(net, seed=run)
        result = simulator.run(
            initial_marking=m0,
            t_max=T_MAX,
            max_steps=MAX_STEPS,
            record_every_step=True,
        )

        times = result.times
        x1 = result.markings[:, x1_index]
        x2 = result.markings[:, x2_index]

        x1_mean = time_weighted_mean(times, x1)
        x2_mean = time_weighted_mean(times, x2)
        x1_slope = initial_slope(times, x1)
        x2_slope = initial_slope(times, x2)
        choice = classify_choice(x1_mean, x2_mean)

        summaries.append(
            {
                "run": run + 1,
                "x1_mean": x1_mean,
                "x2_mean": x2_mean,
                "x1_slope": x1_slope,
                "x2_slope": x2_slope,
                "choice": choice,
                "events": len(result.fired_transitions),
            }
        )

    print(f"Xist-high threshold: {XIST_HIGH_THRESHOLD:.2f}")
    print(f"Initial-slope window: 0-{INITIAL_SLOPE_WINDOW:g} h")
    print("Choice is based on each run's full time-weighted mean.\n")
    print(
        f"{'run':>4}  {'mean x1':>9}  {'mean x2':>9}  "
        f"{'slope x1':>10}  {'slope x2':>10}  "
        f"{'choice':>7}  {'events':>8}"
    )
    for summary in summaries:
        print(
            f"{summary['run']:4d}  {summary['x1_mean']:9.2f}  "
            f"{summary['x2_mean']:9.2f}  {summary['x1_slope']:10.3f}  "
            f"{summary['x2_slope']:10.3f}  {summary['choice']:>7}  "
            f"{summary['events']:8d}"
        )

    choices = ("x1", "x2", "both", "neither")
    choice_counts = {choice: sum(s["choice"] == choice for s in summaries) for choice in choices}
    print("\nChoice counts:")
    for choice in choices:
        print(f"  {choice:7s}: {choice_counts[choice]}")

    runs = np.array([s["run"] for s in summaries])
    width = 0.4
    fig, axes = plt.subplots(3, 1, figsize=(14, 12), constrained_layout=True)

    axes[0].bar(runs - width / 2, [s["x1_mean"] for s in summaries], width, label="x1")
    axes[0].bar(runs + width / 2, [s["x2_mean"] for s in summaries], width, label="x2")
    axes[0].set_title("Time-weighted mean Xist level per run")
    axes[0].set_xlabel("run")
    axes[0].set_ylabel("mean molecule count")
    axes[0].legend()
    axes[0].grid(axis="y", alpha=0.3)

    axes[1].scatter(runs, [s["x1_slope"] for s in summaries], label="x1", s=25)
    axes[1].scatter(runs, [s["x2_slope"] for s in summaries], label="x2", s=25, marker="x")
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_title(f"Initial slopes (first {INITIAL_SLOPE_WINDOW:g} h)")
    axes[1].set_xlabel("run")
    axes[1].set_ylabel("molecules / hour")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    axes[2].bar(choices, [choice_counts[choice] for choice in choices])
    axes[2].set_title(
        "Allele choice from full-run time-weighted means "
        f"(Xist-high threshold = {XIST_HIGH_THRESHOLD:.1f})"
    )
    axes[2].set_xlabel("classification")
    axes[2].set_ylabel("number of runs")
    axes[2].grid(axis="y", alpha=0.3)

    plt.show()
