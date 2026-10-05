"""Reusable pairwise differential-expression analysis with PyDESeq2.

Input format: first CSV column contains feature IDs; remaining columns contain
raw integer counts, one sample per column. Replicate groups are inferred by
removing a trailing ``_<number>`` from each sample name.
"""

from __future__ import annotations

import io
import re
import zipfile
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


@dataclass(frozen=True)
class AnalysisConfig:
    feature_label: str
    min_count: int = 10
    min_samples: int | None = 2
    alpha: float = 0.05
    n_cpus: int = 1


def read_count_matrix(source: str | Path | io.BytesIO) -> pd.DataFrame:
    """Read and validate a feature-by-sample raw-count CSV."""
    frame = pd.read_csv(source, index_col=0)
    if frame.empty or frame.shape[1] < 2:
        raise ValueError("The file must contain features and at least two sample columns.")
    if frame.index.has_duplicates:
        duplicates = frame.index[frame.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Feature IDs must be unique. Duplicates include: {duplicates}")

    numeric = frame.apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        raise ValueError("All count cells must be numeric and non-missing.")
    if (numeric < 0).any().any():
        raise ValueError("Counts cannot be negative.")
    if not np.allclose(numeric.to_numpy(), np.rint(numeric.to_numpy())):
        raise ValueError("DESeq2 requires raw integer counts, not normalized values.")
    return numeric.round().astype(int)


def infer_groups(sample_names: Iterable[str]) -> pd.Series:
    """Infer group labels by stripping a trailing underscore plus digits."""
    names = [str(name) for name in sample_names]
    groups = [re.sub(r"_[0-9]+$", "", name) for name in names]
    metadata = pd.Series(groups, index=names, name="group", dtype="object")
    counts = metadata.value_counts()
    if len(counts) < 2:
        raise ValueError("At least two groups are required after parsing sample names.")
    if (counts < 2).any():
        weak = ", ".join(f"{g} ({n})" for g, n in counts[counts < 2].items())
        raise ValueError(f"Each group needs at least two replicates. Check: {weak}")
    return metadata


def available_comparisons(groups: Iterable[str]) -> list[tuple[str, str]]:
    return list(combinations(sorted(set(groups)), 2))


def fit_dataset(
    counts: pd.DataFrame,
    metadata: pd.Series,
    config: AnalysisConfig,
) -> tuple[DeseqDataSet, int]:
    """Filter low counts and fit the DESeq2-style model once."""
    min_samples = config.min_samples
    if min_samples is None:
        min_samples = int(metadata.value_counts().min())

    keep = (counts >= config.min_count).sum(axis=1) >= min_samples
    filtered = counts.loc[keep]
    if filtered.empty:
        raise ValueError("No features remain after low-count filtering.")

    # PyDESeq2 expects samples in rows and features in columns.
    sample_by_feature = filtered.T
    sample_metadata = pd.DataFrame({"group": metadata.loc[sample_by_feature.index]})
    dds = DeseqDataSet(
        counts=sample_by_feature,
        metadata=sample_metadata,
        design="~group",
        refit_cooks=True,
        n_cpus=config.n_cpus,
        quiet=True,
    )
    dds.deseq2()
    return dds, len(filtered)


def run_pairwise(
    dds: DeseqDataSet,
    pairs: Iterable[tuple[str, str]],
    config: AnalysisConfig,
    progress: Callable[[int, int, str], None] | None = None,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Test selected group pairs; positive LFC means group1 > group2."""
    pairs = list(pairs)
    result_tables: dict[str, pd.DataFrame] = {}
    summary_rows: list[dict[str, object]] = []

    for number, (group1, group2) in enumerate(pairs, start=1):
        name = f"{group1}_vs_{group2}"
        stats = DeseqStats(
            dds,
            contrast=["group", group1, group2],
            alpha=config.alpha,
            n_cpus=config.n_cpus,
            quiet=True,
        )
        stats.summary()

        result = stats.results_df.copy()
        result.index.name = config.feature_label
        result = result.reset_index().sort_values("padj", na_position="last")
        result_tables[name] = result
        significant = int((result["padj"] < config.alpha).fillna(False).sum())
        summary_rows.append(
            {
                "comparison": name,
                "group1": group1,
                "group2": group2,
                "n_tested": len(result),
                f"n_sig_padj{config.alpha:g}": significant,
            }
        )
        if progress:
            progress(number, len(pairs), name)

    sig_col = f"n_sig_padj{config.alpha:g}"
    summary = pd.DataFrame(summary_rows).sort_values(sig_col, ascending=False)
    return result_tables, summary


def save_results(
    output_dir: str | Path,
    result_tables: dict[str, pd.DataFrame],
    summary: pd.DataFrame,
) -> None:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for name, table in result_tables.items():
        table.to_csv(output / f"{name}.txt", sep="\t", index=False)
    summary.to_csv(output / "_summary_all_comparisons.txt", sep="\t", index=False)


def results_zip(
    result_tables: dict[str, pd.DataFrame], summary: pd.DataFrame
) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("_summary_all_comparisons.txt", summary.to_csv(sep="\t", index=False))
        for name, table in result_tables.items():
            archive.writestr(
                f"{name}.txt", table.to_csv(sep="\t", index=False)
            )
    return buffer.getvalue()
