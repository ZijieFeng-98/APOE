#!/usr/bin/env python3
"""High-level clinical analysis workflow for the CGGA APOE cohort.

This module turns the previously sketched analysis outline into an executable
script.  It provides small, composable helpers for summarising the CGGA
clinical table, stratifying by key biomarkers, generating survival plots, and
combining the information with an APOE genotype manifest.  Analysts can run the
complete workflow with a single command or import specific helpers into a
notebook for bespoke exploration.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:  # pragma: no cover - optional heavy imports for type checking
    from lifelines import CoxPHFitter
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

__all__ = [
    "assign_risk_category",
    "create_clinical_summary",
    "stratified_analysis",
    "perform_km_analysis",
    "cox_regression_analysis",
    "analyze_therapy_combinations",
    "subgroup_survival_analysis",
    "merge_apoe_clinical",
    "apoe_survival_analysis",
    "perform_logrank_tests",
    "create_clinical_figure_panel",
    "run_complete_analysis",
    "main",
]


RISK_CATEGORY_MAP = {
    "ε2/ε2": "Protective",
    "ε2/ε3": "Protective",
    "ε2/ε4": "Baseline",
    "ε3/ε3": "Baseline",
    "ε3/ε4": "Increased",
    "ε4/ε4": "High",
}

DEFAULT_APOE_GENOTYPES = {
    "HRR024686": "ε2/ε3",
    "HRR024687": "ε2/ε3",
    "HRR024688": "ε2/ε3",
    "HRR024689": "ε3/ε3",
    "HRR024690": "ε2/ε3",
    "HRR024691": "ε2/ε3",
    "HRR024692": "ε3/ε3",
    "HRR024693": "ε2/ε3",
    "HRR024694": "ε3/ε3",
    "HRR024695": "ε2/ε3",
    "HRR024696": "ε3/ε3",
    "HRR024698": "ε2/ε3",
    "HRR024699": "ε3/ε3",
    "HRR024700": "ε3/ε4",
    "HRR024701": "ε3/ε3",
    "HRR024702": "ε3/ε3",
    "HRR024703": "ε3/ε4",
    "HRR024704": "ε2/ε3",
    "HRR024705": "ε3/ε3",
    "HRR024706": "ε3/ε3",
    "HRR024707": "ε2/ε3",
    "HRR024708": "ε3/ε3",
    "HRR024709": "ε3/ε4",
    "HRR024710": "ε3/ε3",
}


@dataclass(frozen=True)
class AnalysisConfig:
    clinical_path: Path
    output_dir: Path
    apoe_manifest: Optional[Path] = None
    id_mapping: Optional[Path] = None
    time_column: str = "OS"
    event_column: str = "Censor (alive=0; dead=1)"
    group_column: str = "Grade"


def assign_risk_category(genotype: str) -> str:
    """Map an APOE genotype string to a risk category."""

    return RISK_CATEGORY_MAP.get(genotype, "Unknown")


def _require_columns(frame: pd.DataFrame, columns: Sequence[str]) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"DataFrame is missing required columns: {', '.join(missing)}")


def load_apoe_genotypes(path: Optional[Path]) -> pd.DataFrame:
    """Load APOE genotype information from a TSV/CSV file or the defaults."""

    if path is None:
        data = pd.DataFrame(
            sorted(DEFAULT_APOE_GENOTYPES.items()), columns=["Patient_ID", "APOE_genotype"]
        )
    else:
        if not path.exists():
            raise FileNotFoundError(f"APOE genotype manifest not found: {path}")
        if path.suffix.lower() in {".tsv", ".txt"}:
            data = pd.read_csv(path, sep="\t")
        else:
            data = pd.read_csv(path)
    _require_columns(data, ["Patient_ID", "APOE_genotype"])
    data = data.copy()
    data["Risk_Category"] = data["APOE_genotype"].map(assign_risk_category)
    return data


def load_clinical_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Clinical table not found: {path}")
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def create_clinical_summary(df: pd.DataFrame, *, time_column: str, event_column: str) -> pd.DataFrame:
    """Generate cohort-level summary statistics."""

    _require_columns(
        df,
        [
            "Age",
            "Gender",
            time_column,
            event_column,
            "Radio_status (treated=1;un-treated=0)",
            "Chemo_status (TMZ treated=1;un-treated=0)",
            "IDH_mut_status",
            "MGMTp_methylation_status",
        ],
    )

    summary = {
        "Total_Patients": len(df),
        "Age_Mean": df["Age"].mean(),
        "Age_SD": df["Age"].std(),
        "Age_Median": df["Age"].median(),
        "Male_N": (df["Gender"].str.lower() == "male").sum(),
        "Female_N": (df["Gender"].str.lower() == "female").sum(),
        "Time_Median": df[time_column].median(),
        "Death_Rate": df[event_column].mean() * 100,
        "RT_Rate": df["Radio_status (treated=1;un-treated=0)"].mean() * 100,
        "TMZ_Rate": df["Chemo_status (TMZ treated=1;un-treated=0)"].mean() * 100,
        "IDH_Mutant_Rate": (df["IDH_mut_status"].str.lower() == "mutant").mean() * 100,
        "MGMT_Methylated_Rate": (df["MGMTp_methylation_status"].str.lower() == "methylated").mean() * 100,
    }
    return pd.DataFrame([summary])


def stratified_analysis(
    df: pd.DataFrame,
    stratify_by: str,
    *,
    time_column: str,
    event_column: str,
) -> pd.DataFrame:
    """Produce per-group summaries based on ``stratify_by``."""

    _require_columns(df, [stratify_by])
    summaries: List[pd.DataFrame] = []
    for group_value in sorted(df[stratify_by].dropna().unique()):
        subset = df[df[stratify_by] == group_value]
        if subset.empty:
            continue
        summary = create_clinical_summary(subset, time_column=time_column, event_column=event_column)
        summary[stratify_by] = group_value
        summaries.append(summary)
    if not summaries:
        return pd.DataFrame(columns=[stratify_by])
    return pd.concat(summaries, ignore_index=True)


def _prepare_survival_columns(
    df: pd.DataFrame, *, time_column: str, event_column: str
) -> pd.DataFrame:
    _require_columns(df, [time_column, event_column])
    survival_df = df[[time_column, event_column]].copy()
    survival_df = survival_df.dropna(subset=[time_column, event_column])
    survival_df[event_column] = survival_df[event_column].astype(int)
    return survival_df


def perform_km_analysis(
    df: pd.DataFrame,
    group_col: str,
    *,
    time_column: str,
    event_column: str,
    ax: Optional["Axes"] = None,
) -> "Axes":
    """Plot Kaplan–Meier curves and return the axes used for plotting."""

    from lifelines import KaplanMeierFitter
    import matplotlib.pyplot as plt

    _require_columns(df, [group_col])
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))
    kmf = KaplanMeierFitter()
    for group_value in sorted(df[group_col].dropna().unique()):
        mask = df[group_col] == group_value
        subset = df.loc[mask]
        survival_df = _prepare_survival_columns(subset, time_column=time_column, event_column=event_column)
        if survival_df.empty:
            continue
        kmf.fit(
            survival_df[time_column].astype(float),
            event_observed=survival_df[event_column],
            label=str(group_value),
        )
        kmf.plot_survival_function(ax=ax)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Survival Probability")
    ax.set_title(f"Overall Survival by {group_col}")
    ax.legend(loc="best")
    return ax


def cox_regression_analysis(
    df: pd.DataFrame,
    *,
    time_column: str,
    event_column: str,
) -> "CoxPHFitter":
    """Fit a Cox proportional hazards model using common covariates."""

    from lifelines import CoxPHFitter

    required = [
        "Age",
        "Gender",
        "IDH_mut_status",
        "MGMTp_methylation_status",
        "Radio_status (treated=1;un-treated=0)",
        "Chemo_status (TMZ treated=1;un-treated=0)",
    ]
    _require_columns(df, required + [time_column, event_column])
    cox_data = df.copy()
    cox_data["Gender_Male"] = (cox_data["Gender"].str.lower() == "male").astype(int)
    cox_data["IDH_Mutant"] = (cox_data["IDH_mut_status"].str.lower() == "mutant").astype(int)
    cox_data["MGMT_Methylated"] = (cox_data["MGMTp_methylation_status"].str.lower() == "methylated").astype(int)
    cox_data["RT_Treated"] = cox_data["Radio_status (treated=1;un-treated=0)"]
    cox_data["TMZ_Treated"] = cox_data["Chemo_status (TMZ treated=1;un-treated=0)"]

    model_columns = [
        "Age",
        "Gender_Male",
        "IDH_Mutant",
        "MGMT_Methylated",
        "RT_Treated",
        "TMZ_Treated",
    ]
    ready = cox_data[model_columns + [time_column, event_column]].dropna()
    cph = CoxPHFitter()
    cph.fit(ready, duration_col=time_column, event_col=event_column)
    return cph


def analyze_therapy_combinations(
    df: pd.DataFrame,
    *,
    time_column: str,
    event_column: str,
) -> pd.DataFrame:
    """Summarise outcomes by therapy combinations."""

    _require_columns(
        df,
        [
            "Radio_status (treated=1;un-treated=0)",
            "Chemo_status (TMZ treated=1;un-treated=0)",
            time_column,
            event_column,
        ],
    )
    therapy_df = df.copy()
    therapy_df["Therapy_Group"] = "No_Treatment"
    therapy_df.loc[
        (therapy_df["Radio_status (treated=1;un-treated=0)"] == 1)
        & (therapy_df["Chemo_status (TMZ treated=1;un-treated=0)"] == 0),
        "Therapy_Group",
    ] = "RT_Only"
    therapy_df.loc[
        (therapy_df["Radio_status (treated=1;un-treated=0)"] == 0)
        & (therapy_df["Chemo_status (TMZ treated=1;un-treated=0)"] == 1),
        "Therapy_Group",
    ] = "TMZ_Only"
    therapy_df.loc[
        (therapy_df["Radio_status (treated=1;un-treated=0)"] == 1)
        & (therapy_df["Chemo_status (TMZ treated=1;un-treated=0)"] == 1),
        "Therapy_Group",
    ] = "RT_TMZ"

    grouped = therapy_df.groupby("Therapy_Group")
    summary = grouped.agg(
        {
            time_column: ["median", "mean"],
            event_column: "mean",
            "Age": "mean",
            "CGGA_ID": "count",
        }
    )
    summary.columns = ["_".join(filter(None, map(str, col))).strip("_") for col in summary.columns.values]
    return summary.reset_index().rename(columns={f"{event_column}_mean": "Death_Rate"})


def subgroup_survival_analysis(
    df: pd.DataFrame,
    subgroup_col: str,
    subgroup_value: str,
    *,
    time_column: str,
    event_column: str,
) -> Optional["Axes"]:
    """Generate a KM plot for a specific subgroup."""

    _require_columns(df, [subgroup_col])
    subset = df[df[subgroup_col] == subgroup_value]
    if subset.empty:
        return None
    ax = perform_km_analysis(
        subset,
        "Therapy_Group",
        time_column=time_column,
        event_column=event_column,
    )
    ax.set_title(f"Survival in {subgroup_col} = {subgroup_value}")
    return ax


def merge_apoe_clinical(
    apoe_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    *,
    id_mapping: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Merge APOE genotype information with the CGGA clinical table."""

    apoe_df = apoe_df.copy()
    _require_columns(apoe_df, ["Patient_ID", "APOE_genotype", "Risk_Category"])
    if id_mapping is not None:
        _require_columns(id_mapping, ["Patient_ID", "CGGA_ID"])
        merged = id_mapping.merge(apoe_df, on="Patient_ID", how="inner")
        return merged.merge(clinical_df, on="CGGA_ID", how="inner")
    return apoe_df.merge(clinical_df, left_on="Patient_ID", right_on="CGGA_ID", how="inner")


def apoe_survival_analysis(
    merged_df: pd.DataFrame,
    *,
    time_column: str,
    event_column: str,
) -> List["Axes"]:
    """Return Kaplan–Meier plots for APOE genotype and risk category."""

    axes: List["Axes"] = []
    for column in ["APOE_genotype", "Risk_Category"]:
        if column in merged_df.columns:
            axes.append(
                perform_km_analysis(
                    merged_df,
                    column,
                    time_column=time_column,
                    event_column=event_column,
                )
            )
    return axes


def perform_logrank_tests(
    df: pd.DataFrame,
    group_col: str,
    *,
    time_column: str,
    event_column: str,
) -> pd.DataFrame:
    """Run pairwise log-rank tests between groups."""

    from lifelines.statistics import logrank_test

    _require_columns(df, [group_col])
    groups = sorted(df[group_col].dropna().unique())
    results = []
    for i, g1 in enumerate(groups):
        for g2 in groups[i + 1 :]:
            mask1 = df[group_col] == g1
            mask2 = df[group_col] == g2
            surv1 = _prepare_survival_columns(df.loc[mask1], time_column=time_column, event_column=event_column)
            surv2 = _prepare_survival_columns(df.loc[mask2], time_column=time_column, event_column=event_column)
            if surv1.empty or surv2.empty:
                continue
            test = logrank_test(
                surv1[time_column].astype(float),
                surv2[time_column].astype(float),
                event_observed_A=surv1[event_column],
                event_observed_B=surv2[event_column],
            )
            results.append({"Group1": g1, "Group2": g2, "Test_Statistic": test.test_statistic, "P_Value": test.p_value})
    return pd.DataFrame(results)


def create_clinical_figure_panel(df: pd.DataFrame, *, time_column: str) -> "Figure":
    """Create a multi-panel overview figure of key cohort attributes."""

    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    sns.histplot(df, x="Age", bins=20, ax=axes[0, 0])
    axes[0, 0].set_title("Age Distribution")

    gender_counts = df["Gender"].value_counts(dropna=False)
    axes[0, 1].pie(gender_counts.values, labels=gender_counts.index, autopct="%1.1f%%")
    axes[0, 1].set_title("Gender Distribution")

    therapy_totals = df[[
        "Radio_status (treated=1;un-treated=0)",
        "Chemo_status (TMZ treated=1;un-treated=0)",
    ]].sum()
    axes[0, 2].bar(therapy_totals.index, therapy_totals.values)
    axes[0, 2].set_title("Therapy Distribution")

    sns.countplot(data=df, x="IDH_mut_status", ax=axes[1, 0])
    axes[1, 0].set_title("IDH Mutation Status")

    sns.countplot(data=df, x="MGMTp_methylation_status", ax=axes[1, 1])
    axes[1, 1].set_title("MGMT Methylation Status")
    axes[1, 1].tick_params(axis="x", rotation=45)

    sns.histplot(df, x=time_column, bins=30, ax=axes[1, 2])
    axes[1, 2].set_title("Survival Time Distribution")
    axes[1, 2].set_xlabel(time_column)

    fig.tight_layout()
    return fig


def _load_optional_mapping(path: Optional[Path]) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    if not path.exists():
        raise FileNotFoundError(f"ID mapping file not found: {path}")
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def run_complete_analysis(config: AnalysisConfig) -> None:
    """Execute the full analysis pipeline and persist results."""

    import matplotlib.pyplot as plt

    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    clinical_df = load_clinical_table(config.clinical_path)
    clinical_df = clinical_df.copy()

    apoe_df = load_apoe_genotypes(config.apoe_manifest)
    id_mapping_df = _load_optional_mapping(config.id_mapping)

    # Cohort-level summaries
    overall_summary = create_clinical_summary(
        clinical_df, time_column=config.time_column, event_column=config.event_column
    )
    stratified_summary = stratified_analysis(
        clinical_df,
        config.group_column,
        time_column=config.time_column,
        event_column=config.event_column,
    )

    km_axes = perform_km_analysis(
        clinical_df,
        config.group_column,
        time_column=config.time_column,
        event_column=config.event_column,
    )
    km_figure = km_axes.get_figure()
    km_path = output_dir / f"survival_by_{config.group_column}.png"
    km_figure.savefig(km_path, dpi=300, bbox_inches="tight")
    plt.close(km_figure)

    cph = cox_regression_analysis(
        clinical_df, time_column=config.time_column, event_column=config.event_column
    )

    therapy_summary = analyze_therapy_combinations(
        clinical_df, time_column=config.time_column, event_column=config.event_column
    )

    logrank_results = perform_logrank_tests(
        clinical_df,
        config.group_column,
        time_column=config.time_column,
        event_column=config.event_column,
    )

    panel_fig = create_clinical_figure_panel(clinical_df, time_column=config.time_column)
    panel_path = output_dir / "clinical_overview_panel.png"
    panel_fig.savefig(panel_path, dpi=300, bbox_inches="tight")
    plt.close(panel_fig)

    merged_df = merge_apoe_clinical(apoe_df, clinical_df, id_mapping=id_mapping_df)
    apoe_axes = apoe_survival_analysis(
        merged_df,
        time_column=config.time_column,
        event_column=config.event_column,
    )
    apoe_paths = []
    for ax, label in zip(
        apoe_axes,
        ["apoe_genotype", "apoe_risk_category"],
    ):
        fig = ax.get_figure()
        path = output_dir / f"survival_by_{label}.png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        apoe_paths.append(path)
        plt.close(fig)

    # Persist tabular outputs
    results_book = output_dir / "clinical_analysis_results.xlsx"
    with pd.ExcelWriter(results_book) as writer:
        overall_summary.to_excel(writer, sheet_name="Overall_Summary", index=False)
        stratified_summary.to_excel(writer, sheet_name=f"By_{config.group_column}", index=False)
        therapy_summary.to_excel(writer, sheet_name="Therapy_Analysis", index=False)
        logrank_results.to_excel(writer, sheet_name="LogRank_Tests", index=False)
        merged_df.to_excel(writer, sheet_name="Merged_APOE", index=False)
        cph.summary.to_excel(writer, sheet_name="Cox_Model")

    manifest = {
        "overall_summary": str(results_book),
        "survival_by_group": str(km_path),
        "overview_panel": str(panel_path),
        "cox_model": str(results_book),
        "therapy_summary": str(results_book),
        "logrank_tests": str(results_book),
        "merged_cohort": str(results_book),
        "apoe_survival_plots": [str(path) for path in apoe_paths],
    }
    manifest_path = output_dir / "clinical_analysis_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True, type=Path, help="CGGA clinical data table (TSV/CSV).")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("clinical_analysis_outputs"),
        help="Directory to store generated reports.",
    )
    parser.add_argument(
        "--apoe-genotypes",
        type=Path,
        help="Optional TSV/CSV file with Patient_ID and APOE_genotype columns.",
    )
    parser.add_argument(
        "--id-mapping",
        type=Path,
        help="Optional mapping file with Patient_ID and CGGA_ID columns.",
    )
    parser.add_argument(
        "--time-column",
        default="OS",
        help="Clinical time-to-event column (default: OS).",
    )
    parser.add_argument(
        "--event-column",
        default="Censor (alive=0; dead=1)",
        help="Event indicator column (default: Censor (alive=0; dead=1)).",
    )
    parser.add_argument(
        "--group-column",
        default="Grade",
        help="Column used for default stratified analyses (default: Grade).",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    config = AnalysisConfig(
        clinical_path=args.clinical.resolve(),
        output_dir=args.output_dir.resolve(),
        apoe_manifest=args.apoe_genotypes.resolve() if args.apoe_genotypes else None,
        id_mapping=args.id_mapping.resolve() if args.id_mapping else None,
        time_column=args.time_column,
        event_column=args.event_column,
        group_column=args.group_column,
    )
    run_complete_analysis(config)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry-point
    raise SystemExit(main())

