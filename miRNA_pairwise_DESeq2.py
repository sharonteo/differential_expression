"""Python replacement for miRNA_pairwise_DESeq2.R."""

import argparse

from deseq2_pairwise import (
    AnalysisConfig,
    available_comparisons,
    fit_dataset,
    infer_groups,
    read_count_matrix,
    run_pairwise,
    save_results,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("counts_csv")
    parser.add_argument("--output", default="DE_results_miRNA")
    parser.add_argument("--cpus", type=int, default=1)
    args = parser.parse_args()

    counts = read_count_matrix(args.counts_csv)
    metadata = infer_groups(counts.columns)
    config = AnalysisConfig(
        feature_label="miRNA",
        min_samples=2,
        n_cpus=args.cpus,
    )
    dds, retained = fit_dataset(counts, metadata, config)
    print(f"Detected groups:\n{metadata.value_counts().sort_index().to_string()}")
    print(f"Kept {retained} of {len(counts)} miRNAs after filtering")
    tables, summary = run_pairwise(dds, available_comparisons(metadata), config)
    save_results(args.output, tables, summary)
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
