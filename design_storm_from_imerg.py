"""
Design storms from IMERG half-hourly events + daily frequency analysis.

1. Reads the half-hourly event CSVs (Giovanni area-average, mm/hr).
2. For each event finds the maximum depth for durations 0.5-24 h and divides by the event's
   maximum 24-h depth  ->  depth-duration RATIOS.  The mean ratio curve is the storm shape.
3. 24-h design depth = Gumbel daily depth (from rainfall_analysis.py) x f24, where f24 is the mean ratio
   (max 24 h)/(max UTC day) measured on the same events (replaces the textbook 1.13).
4. Alternating-block hyetographs (30-min blocks) for each return period -> SWMM TIMESERIES (mm/h).

Usage:  python design_storm_from_imerg.py --events DIR --depths rain_output/return_period_depths.csv --out design_storms
Requires numpy, pandas, scipy, matplotlib
"""
import argparse, glob, os, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DUR_H = np.array([0.5, 1, 2, 3, 6, 12, 24])
STEPS = (DUR_H * 2).astype(int)          # number of 30-min steps


def load_event(path):
    df = pd.read_csv(path, skiprows=8)
    df.columns = ["time", "rate"]
    df = df.dropna()
    df["time"] = pd.to_datetime(df["time"])
    df.loc[df["rate"] < -999, "rate"] = np.nan
    if df["rate"].isna().any():
        raise ValueError("fill values in " + path)
    return df.set_index("time")["rate"] * 0.5            # mm/h -> mm per 30 min


def event_stats(depth):
    d = depth.values
    maxd = np.array([np.convolve(d, np.ones(k), "valid").max() for k in STEPS])
    utc_day = depth.groupby(depth.index.date).sum()
    return maxd, utc_day.max(), utc_day.idxmax()


def alternating_block(cum_fn, total_h=24, dt_h=0.5):
    n = int(round(total_h / dt_h))
    cum = np.array([cum_fn(dt_h * i) for i in range(n + 1)])
    blocks = np.diff(cum)
    order = np.argsort(-blocks)
    out = np.zeros(n); peak = n // 2
    out[peak] = blocks[order[0]]; left, right = peak - 1, peak + 1; go_right = True
    for k in order[1:]:
        if (go_right and right < n) or left < 0:
            out[right] = blocks[k]; right += 1
        else:
            out[left] = blocks[k]; left -= 1
        go_right = not go_right
    return out                                              # mm per block


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True); ap.add_argument("--depths", required=True)
    ap.add_argument("--out", default="design_storms"); ap.add_argument("--high-b", type=float, default=0.6, help="exponent b of the high-intensity scenario D(d)=D24*(d/24)^(1-b)"); ap.add_argument("--daily", help="optional daily CSV (date,rain_mm) for a consistency check")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.events, "*GPM_3IMERGHH*.csv")))
    rows, ratios, curves = [], [], {}
    daily = pd.read_csv(a.daily, parse_dates=["date"]).set_index("date")["rain_mm"] if a.daily else None
    for f in files:
        dep = load_event(f); maxd, dmax, dday = event_stats(dep)
        name = re.search(r"(\d{8})-(\d{8})", f).group(1)
        r = maxd / maxd[-1]; ratios.append(r); curves[name] = dep
        row = {"event_window_start": name, "max_UTC_day_mm": dmax, "UTC_day": str(dday), "max_24h_mm": maxd[-1],
               "f24": maxd[-1] / dmax, "peak_30min_mm": maxd[0], "peak_1h_mm": maxd[1], "peak_3h_mm": maxd[3]}
        if daily is not None: row["daily_file_mm"] = daily.get(pd.Timestamp(dday), np.nan)
        rows.append(row)
    ev = pd.DataFrame(rows).round(2); ev.to_csv(os.path.join(a.out, "event_summary.csv"), index=False)
    R = np.array(ratios); mean_r = R.mean(0)
    ddf = pd.DataFrame({"duration_h": DUR_H, "mean_ratio": mean_r, "min": R.min(0), "max": R.max(0)}).round(3)
    ddf.to_csv(os.path.join(a.out, "depth_duration_ratios.csv"), index=False)
    f24 = ev["f24"].mean()
    print(ev.to_string(index=False)); print("\nmean f24 = %.3f\n" % f24); print(ddf.to_string(index=False))
    # smooth cumulative curve: log-log linear interpolation through the mean ratios (D(0)=0)
    ld, lr = np.log(DUR_H), np.log(mean_r)
    cum_ratio = lambda t: 0.0 if t <= 0 else float(np.exp(np.interp(np.log(t), ld, lr)))
    dep = pd.read_csv(a.depths); lines = []; storms = {}
    for _, row in dep.iterrows():
        T = int(row["T_years"])
        if T not in (2, 5, 10, 25, 50): continue
        d24 = row["Gumbel_mm"] * f24
        blocks = alternating_block(lambda t: d24 * cum_ratio(t))
        inten = blocks / 0.5; storms[T] = (d24, inten)
        for i, v in enumerate(inten): lines.append("T%-4d %02d:%02d  %.3f\n" % (T, (i * 30) // 60, (i * 30) % 60, v))
        lines.append("T%-4d 24:00  0.000\n" % T)
        print("T=%d y: 24-h depth %.1f mm | peak 30-min %.1f mm/h | peak 1-h %.1f mm" % (T, d24, inten.max(), np.convolve(blocks, [1, 1], 'valid').max()))
    # high-intensity SENSITIVITY scenario (assumption, not data): power-law depth-duration curve
    for T in (2, 10, 50):
        d24 = storms[T][0]; blocks = alternating_block(lambda t: d24 * (t / 24.0) ** (1 - a.high_b)); inten = blocks / 0.5
        storms["H%d" % T] = (d24, inten)
        for i, v in enumerate(inten): lines.append("H%-4d %02d:%02d  %.3f\n" % (T, (i * 30) // 60, (i * 30) % 60, v))
        lines.append("H%-4d 24:00  0.000\n" % T)
        print("H%d (b=%.1f): 24-h %.1f mm | peak 30-min %.1f mm/h | peak 1-h %.1f mm" % (T, a.high_b, d24, inten.max(), np.convolve(blocks, [1, 1], 'valid').max()))
    open(os.path.join(a.out, "design_storms_swmm.txt"), "w").write(";; Design storms, INTENSITY mm/h, 30-min steps, 24 h, peak at 12:00\n" + "".join(lines))
    # plots
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    for r in R: ax[0].plot(DUR_H, r, color="grey", alpha=0.5, lw=1)
    ax[0].plot(DUR_H, mean_r, "k-o", lw=2, label="mean of 8 events"); ax[0].set_xscale("log"); ax[0].set_xlabel("Duration (h)")
    ax[0].set_ylabel("Max depth / max 24-h depth"); ax[0].grid(alpha=0.3, which="both"); ax[0].legend(); ax[0].set_title("Depth-duration ratios (IMERG half-hourly)")
    for T, (d24, inten) in storms.items():
        if isinstance(T, int): ax[1].step(np.arange(48) * 0.5, inten, where="post", label="T=%d y (%.0f mm)" % (T, d24))
    ax[1].step(np.arange(48) * 0.5, storms["H10"][1], where="post", color="k", ls="--", label="H10: high-intensity sensitivity")
    ax[1].set_xlabel("Time (h)"); ax[1].set_ylabel("Intensity (mm/h)"); ax[1].grid(alpha=0.3); ax[1].legend(); ax[1].set_title("Design storms (24 h, 30-min blocks)")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "design_storms.png"), dpi=170)


if __name__ == "__main__":
    main()
