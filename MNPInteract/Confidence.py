"""Confidence.py — Stage 5: Final high-confidence interactor generation."""

import pandas as pd
from pathlib import Path

DEFAULT_GO_TERMS = [
    "plasmodesma",
    "membrane",
    "plasma membrane",
    "cell wall",
    "apoplast",
    "extracellular region",
]


class Confidence(object):

    def __init__(self, config):
        self.config = config
        out = config.output_dir

        self.annotated_csv = out / "Prioritized candidate list_with_PFAM_GO_annotations.csv"
        self.topology_csv  = out / "topology_contact_report.csv"

        self.full_output       = out / "Prioritized_candidates_with_high_confidence.csv"
        self.summary_output    = out / "high_confidence_summary.csv"
        self.true_output       = out / "High_Confidence_Interactors_only.csv"
        self.final_true_output = out / "Final_High_Confidence_Interactors_only.csv"

        self.iptm_cutoff  = config.iptm_cutoff
        self.dockq_cutoff = config.dockq_cutoff
        self.go_terms     = DEFAULT_GO_TERMS

    def run(self):
        print("\n" + "="*60, flush=True)
        print("STAGE 5 — Final high-confidence interactor generation", flush=True)
        print("="*60, flush=True)

        self._check_inputs()

        df_prioritized = pd.read_csv(self.annotated_csv)
        df_topology    = pd.read_csv(self.topology_csv)

        df_prioritized.columns = df_prioritized.columns.str.strip()
        df_topology.columns    = df_topology.columns.str.strip()
        df_prioritized["jobs"] = df_prioritized["jobs"].astype(str).str.strip()
        df_topology["folder"]  = df_topology["folder"].astype(str).str.strip()

        # Normalise mpDockQ column name
        if "mpDockQ/pDockQ" in df_prioritized.columns and "mpDockQ.pDockQ" not in df_prioritized.columns:
            df_prioritized = df_prioritized.rename(columns={"mpDockQ/pDockQ": "mpDockQ.pDockQ"})

        df = df_prioritized.merge(
            df_topology[["folder", "interaction_verdict"]],
            left_on="jobs", right_on="folder", how="left"
        )
        if "folder" in df.columns:
            df.drop(columns=["folder"], inplace=True)

        df = self._reorder_verdict(df)

        df["go_cc_term"]          = df["go_cc_term"].fillna("").astype(str).str.strip()
        df["interaction_verdict"] = df["interaction_verdict"].fillna("").astype(str).str.strip()
        df["iptm"]                = pd.to_numeric(df["iptm"], errors="coerce")
        df["mpDockQ.pDockQ"]      = pd.to_numeric(df["mpDockQ.pDockQ"], errors="coerce")

        df["High Confidence Interactors"] = df.apply(
            lambda row: "TRUE"
            if self._go_matches(row["go_cc_term"])
            and row["interaction_verdict"].lower() == "likely"
            else "FALSE",
            axis=1
        )

        df.to_csv(self.full_output, index=False)
        print("[Confidence] Full dataset   : {}".format(self.full_output), flush=True)

        summary_df = self._make_summary(df)
        summary_df.to_csv(self.summary_output, index=False)
        print("[Confidence] Summary        : {}".format(self.summary_output), flush=True)

        true_df = df[df["High Confidence Interactors"] == "TRUE"].copy()
        true_df.to_csv(self.true_output, index=False)
        print("[Confidence] HC interactors : {}".format(self.true_output), flush=True)

        final_df = df[
            (df["High Confidence Interactors"] == "TRUE") &
            (df["iptm"]           >= self.iptm_cutoff) &
            (df["mpDockQ.pDockQ"] >= self.dockq_cutoff)
        ].copy()
        final_df.to_csv(self.final_true_output, index=False)
        print("[Confidence] Final HC list  : {}".format(self.final_true_output), flush=True)

        print("\n" + "-"*60, flush=True)
        print("RESULTS", flush=True)
        print("-"*60, flush=True)
        print("  Total candidates merged   : {}".format(len(df)), flush=True)
        print("  High Confidence (TRUE)    : {}".format(len(true_df)), flush=True)
        print("  Final HC (+ score filter) : {}".format(len(final_df)), flush=True)
        print("\n  Verdict counts:", flush=True)
        print("  {}".format(summary_df.to_string(index=False)), flush=True)

        return final_df

    def _check_inputs(self):
        missing = []
        if not self.annotated_csv.exists():
            missing.append(str(self.annotated_csv))
        if not self.topology_csv.exists():
            missing.append(str(self.topology_csv))
        if missing:
            raise FileNotFoundError(
                "Required files not found:\n" +
                "\n".join("  {}".format(p) for p in missing) +
                "\nEnsure Stages 1-4 completed successfully."
            )

    def _go_matches(self, text):
        if text == "":
            return True
        text_lower = text.lower()
        return any(term in text_lower for term in self.go_terms)

    @staticmethod
    def _reorder_verdict(df):
        cols = list(df.columns)
        if "interaction_verdict" in cols and "go_cc_term" in cols:
            cols.remove("interaction_verdict")
            idx = cols.index("go_cc_term") + 1
            cols.insert(idx, "interaction_verdict")
            df = df[cols]
        return df

    @staticmethod
    def _make_summary(df):
        summary = (
            df["High Confidence Interactors"]
            .value_counts()
            .rename_axis("High Confidence Interactors")
            .reset_index(name="count")
        )
        base = pd.DataFrame({"High Confidence Interactors": ["TRUE", "FALSE"]})
        summary = base.merge(summary, on="High Confidence Interactors", how="left").fillna(0)
        summary["count"] = summary["count"].astype(int)
        return summary
