import re
import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.ndimage import gaussian_filter1d

def parse_bin_size(size_token, bin_size_map):
    if size_token in bin_size_map:
        return int(bin_size_map[size_token])

    m = re.search(r"(\d+)", str(size_token))
    if not m:
        raise ValueError(f"Could not infer bin size from wildcard: {size_token}")
    return int(m.group(1))


def read_rt_bg(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "End", "RT"],
        dtype={"Chromosome": str}
    )
    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df["RT"] = pd.to_numeric(df["RT"], errors="coerce")
    df = df.dropna(subset=["Chromosome", "Start", "End", "RT"]).copy()
    df = df.sort_values(["Chromosome", "Start"]).reset_index(drop=True)
    return df


def read_gap_file(path):
    if not path:
        return None
    gaps = pd.read_csv(path, sep="\t")
    required = {"Chromosome", "Start", "End"}
    if not required.issubset(gaps.columns):
        raise ValueError(f"Gap file must contain {required}")
    gaps["Chromosome"] = gaps["Chromosome"].astype(str)
    gaps["Start"] = pd.to_numeric(gaps["Start"], errors="coerce")
    gaps["End"] = pd.to_numeric(gaps["End"], errors="coerce")
    gaps = gaps.dropna(subset=["Chromosome", "Start", "End"]).copy()
    return gaps


def read_target_distribution(path):
    if not path:
        raise ValueError("Quantile normalization requested but no target_file provided.")
    target = pd.read_csv(path, sep="\t")
    if "Mean" not in target.columns:
        raise ValueError("Target file must contain a 'Mean' column.")
    arr = pd.to_numeric(target["Mean"], errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        raise ValueError("Target distribution contains no valid numeric values.")
    return arr


def robust_center(values):
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values)
    return values - med


def robust_iqr_scale(values, do_center=True):
    values = np.asarray(values, dtype=float)
    if do_center:
        values = values - np.nanmedian(values)
    q75, q25 = np.nanpercentile(values, [75, 25])
    iqr = q75 - q25
    if not np.isfinite(iqr) or iqr == 0:
        raise ValueError("IQR is zero or non-finite; cannot perform IQR scaling.")
    return values * 1.59 / iqr


def quantile_normalize_to_target(values, target):
    values = np.asarray(values, dtype=float)
    target = np.asarray(target, dtype=float)

    mask = np.isfinite(values)
    out = np.full(values.shape, np.nan, dtype=float)
    x = values[mask]

    if x.size == 0:
        return out

    target_sorted = np.sort(target[np.isfinite(target)])
    if target_sorted.size == 0:
        raise ValueError("Target distribution has no finite values.")

    if target_sorted.size != x.size:
        src_q = np.linspace(0, 1, target_sorted.size)
        dst_q = np.linspace(0, 1, x.size)
        target_rank_values = np.interp(dst_q, src_q, target_sorted)
    else:
        target_rank_values = target_sorted

    order = np.argsort(x, kind="mergesort")
    normalized = np.empty_like(x, dtype=float)
    normalized[order] = target_rank_values
    out[mask] = normalized
    return out


def apply_normalization(df, method, median_center, target=None):
    vals = df["RT"].to_numpy(dtype=float)

    if method == "none":
        out = vals.copy()

    elif method == "median":
        out = robust_center(vals) if median_center else vals.copy()

    elif method == "iqr":
        out = robust_iqr_scale(vals, do_center=median_center)

    elif method == "quantile":
        pre = robust_center(vals) if median_center else vals.copy()
        out = quantile_normalize_to_target(pre, target)

    else:
        raise ValueError(f"Unsupported normalization method: {method}")

    df = df.copy()
    df["RT_norm"] = out
    return df


def loess_smooth_per_chromosome(
    df,
    span_bp=500000,
    exclude_chrY=True,
    clip_min=-8,
    clip_max=8,
):
    df = df.copy()
    df["RT_smooth"] = np.nan

    chrs = df["Chromosome"].unique().tolist()
    if exclude_chrY:
        chrs = [c for c in chrs if c != "chrY"]

    for chrom in chrs:
        sub_idx = df["Chromosome"] == chrom
        sub = df.loc[sub_idx, ["Start", "RT_norm"]].dropna().sort_values("Start")
        if sub.empty:
            continue

        if len(sub) < 3:
            smoothed = sub["RT_norm"].to_numpy(dtype=float)
            df.loc[sub.index, "RT_smooth"] = smoothed
            continue

        x = sub["Start"].to_numpy(dtype=float)
        y = sub["RT_norm"].to_numpy(dtype=float)

        chrom_span = x.max() - x.min()
        if chrom_span <= 0:
            frac = 1.0
        else:
            frac = span_bp / chrom_span

        frac = min(max(frac, 1.0 / len(sub)), 1.0)

        smoothed = lowess(
            endog=y,
            exog=x,
            frac=frac,
            it=0,
            return_sorted=False,
        )

        smoothed[(smoothed <= clip_min) | (smoothed >= clip_max)] = np.nan
        df.loc[sub.index, "RT_smooth"] = smoothed

    if not exclude_chrY and "chrY" in df["Chromosome"].values:
        y_idx = (df["Chromosome"] == "chrY") & (df["RT_smooth"].isna())
        df.loc[y_idx, "RT_smooth"] = df.loc[y_idx, "RT_norm"]

    return df



def gaussian_smooth_per_chromosome(
    df,
    sigma_bp=100000,
    exclude_chrY=True,
    clip_min=-8,
    clip_max=8,
    value_col="RT_norm",
    output_col="RT_smooth",
):
    """
    Gaussian smoothing for regularly binned RT data.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain Chromosome, Start, End, and value_col.
    sigma_bp : int
        Gaussian sigma in base pairs.
    exclude_chrY : bool
        Whether to skip chrY.
    clip_min, clip_max : float
        Values outside this range after smoothing are set to NaN.
    value_col : str
        Column to smooth.
    output_col : str
        Output smoothed column name.
    """
    df = df.copy()
    df[output_col] = np.nan

    chrs = df["Chromosome"].unique().tolist()
    if exclude_chrY:
        chrs = [c for c in chrs if c != "chrY"]

    for chrom in chrs:
        sub_idx = df["Chromosome"] == chrom
        sub = df.loc[sub_idx].sort_values("Start").copy()

        if sub.empty:
            continue

        x = sub[value_col].to_numpy(dtype=float)

        # infer bin size from genomic coordinates
        if len(sub) > 1:
            diffs = np.diff(sub["Start"].to_numpy(dtype=float))
            diffs = diffs[diffs > 0]
            if len(diffs) == 0:
                sigma_bins = 1.0
            else:
                bin_size = np.median(diffs)
                sigma_bins = max(float(sigma_bp) / float(bin_size), 1e-6)
        else:
            sigma_bins = 1.0

        valid = np.isfinite(x)
        if valid.sum() == 0:
            continue

        # NaN-aware Gaussian smoothing:
        # smooth the values with NaNs filled as 0,
        # smooth the valid-mask,
        # divide to reweight properly.
        x_filled = np.where(valid, x, 0.0)
        weights = valid.astype(float)

        smooth_num = gaussian_filter1d(
            x_filled,
            sigma=sigma_bins,
            mode="nearest",
            truncate=4.0,
        )
        smooth_den = gaussian_filter1d(
            weights,
            sigma=sigma_bins,
            mode="nearest",
            truncate=4.0,
        )

        smoothed = np.divide(
            smooth_num,
            smooth_den,
            out=np.full_like(smooth_num, np.nan, dtype=float),
            where=smooth_den > 1e-12,
        )

        smoothed[(smoothed <= clip_min) | (smoothed >= clip_max)] = np.nan
        df.loc[sub.index, output_col] = smoothed

    if not exclude_chrY and "chrY" in df["Chromosome"].values:
        y_idx = (df["Chromosome"] == "chrY") & (df[output_col].isna())
        df.loc[y_idx, output_col] = df.loc[y_idx, value_col]

    return df

def mask_gaps(df, gaps):
    if gaps is None:
        return df

    df = df.copy()
    for _, gap in gaps.iterrows():
        chrom = gap["Chromosome"]
        start = gap["Start"]
        end = gap["End"]
        m = (
            (df["Chromosome"] == chrom)
            & (df["Start"] >= start)
            & (df["Start"] <= end)
        )
        df.loc[m, "RT_smooth"] = np.nan
    return df


def main():
    inp = snakemake.input.rt
    outp = snakemake.output.rt

    norm_method = snakemake.params.norm_method
    median_center = bool(snakemake.params.median_center)
    target_file = snakemake.params.target_file
    gap_file = snakemake.params.gap_file
    exclude_chrY = bool(snakemake.params.exclude_chrY)
    loess_span_bp = int(snakemake.params.loess_span_bp)
    mask_gaps_flag = bool(snakemake.params.mask_gaps)
    clip_min = float(snakemake.params.clip_min)
    clip_max = float(snakemake.params.clip_max)
    bin_size_map = dict(snakemake.params.bin_size_map)
    smoothing = snakemake.params.smooth

    size_token = snakemake.wildcards.sizes
    bin_size = parse_bin_size(size_token, bin_size_map)

    # Currently used only for logging / potential extensions
    span_bins = max(1, int(round(loess_span_bp / bin_size)))
    print(f"[normalize_smooth_rt] input={inp}")
    print(f"[normalize_smooth_rt] normalization={norm_method}")
    print(f"[normalize_smooth_rt] median_center={median_center}")
    print(f"[normalize_smooth_rt] size_token={size_token}, bin_size={bin_size}, span_bp={loess_span_bp}, approx_span_bins={span_bins}")

    df = read_rt_bg(inp)

    target = None
    if norm_method == "quantile":
        target = read_target_distribution(target_file)

    gaps = read_gap_file(gap_file) if gap_file else None

    df = apply_normalization(
        df=df,
        method=norm_method,
        median_center=median_center,
        target=target,
    )
    if smoothing == 'loess':
        df = loess_smooth_per_chromosome(
            df=df,
            span_bp=loess_span_bp,
            exclude_chrY=exclude_chrY,
            clip_min=clip_min,
            clip_max=clip_max,
        )
    if smoothing == 'gaussian':
        df = gaussian_smooth_per_chromosome(
            df=df,
            sigma_bp=loess_span_bp,
            exclude_chrY=exclude_chrY,
            clip_min=clip_min,
            clip_max=clip_max,
            value_col="RT_norm",
            output_col="RT_smooth",
        )
    if mask_gaps_flag:
        df = mask_gaps(df, gaps)

    out = df[["Chromosome", "Start", "End", "RT_smooth"]].dropna(subset=["RT_smooth"]).copy()
    out.to_csv(outp, sep="\t", header=False, index=False, na_rep="NA")


if __name__ == "__main__":
    main()
