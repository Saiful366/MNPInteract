"""Annotate.py — Stage 1: Score filtering and PFAM/GO annotation."""

import time
import requests
import pandas as pd
from pathlib import Path

try:
    from pybiomart import Dataset as BiomartDataset
except ImportError:
    BiomartDataset = None


class Annotate(object):

    def __init__(self, config):
        self.config = config
        out = config.output_dir

        self.input_csv     = config.interpae_csv
        self.cutoff_csv    = out / "Prioritized candidate list.csv"
        self.annotated_csv = out / "Prioritized candidate list_with_PFAM_GO_annotations.csv"
        self.pfam_csv      = out / "pfam_summary_with_ATids_and_annotations.csv"

        self.iptm_cutoff    = config.iptm_cutoff
        self.dockq_cutoff   = config.dockq_cutoff
        self.piscore_cutoff = config.piscore_cutoff

    def run(self):
        if self.annotated_csv.exists():
            print("[Annotate] Skipping — annotated CSV already exists:", flush=True)
            print("  {}".format(self.annotated_csv), flush=True)
            return self.annotated_csv

        print("\n" + "="*60, flush=True)
        print("STAGE 1 — Score filtering and PFAM/GO annotation", flush=True)
        print("="*60, flush=True)

        df = self._load_and_filter()
        df = self._add_gene_id(df)
        df = self._annotate_biomart(df)
        df = self._annotate_interpro(df)
        df = self._finalise_columns(df)

        df.to_csv(self.annotated_csv, index=False)
        print("\n[Annotate] Annotated CSV saved: {}".format(self.annotated_csv), flush=True)
        return self.annotated_csv

    def _load_and_filter(self):
        print("[Annotate] Reading: {}".format(self.input_csv), flush=True)
        df = pd.read_csv(self.input_csv)

        # Normalise mpDockQ column name — AlphaPulldown versions use either
        # 'mpDockQ/pDockQ' (slash) or 'mpDockQ.pDockQ' (dot)
        if "mpDockQ/pDockQ" in df.columns and "mpDockQ.pDockQ" not in df.columns:
            df = df.rename(columns={"mpDockQ/pDockQ": "mpDockQ.pDockQ"})
            print("[Annotate] Renamed column 'mpDockQ/pDockQ' -> 'mpDockQ.pDockQ'", flush=True)

        required = ["jobs", "iptm", "mpDockQ.pDockQ", "pi_score"]
        missing  = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(
                "Input CSV missing required columns: {}\nFound: {}".format(
                    missing, list(df.columns))
            )

        df["iptm"]           = pd.to_numeric(df["iptm"],           errors="coerce")
        df["mpDockQ.pDockQ"] = pd.to_numeric(df["mpDockQ.pDockQ"], errors="coerce")
        df["pi_score"]       = pd.to_numeric(df["pi_score"],       errors="coerce")

        df_filtered = df[
            (df["iptm"]           >= self.iptm_cutoff)  |
            (df["mpDockQ.pDockQ"] >= self.dockq_cutoff) |
            (df["pi_score"]       >  self.piscore_cutoff)
        ].copy()

        df_filtered.to_csv(self.cutoff_csv, index=False)
        print("[Annotate] Cutoff filter: {} -> {} rows".format(len(df), len(df_filtered)), flush=True)
        print("           Saved: {}".format(self.cutoff_csv), flush=True)
        return df_filtered

    def _add_gene_id(self, df):
        df = df.copy()
        df["jobs_new"] = (
            df["jobs"].astype(str)
            .str.replace("PDLP5_and_", "", regex=False)
            .str.strip()
        )
        return df

    def _annotate_biomart(self, df):
        if BiomartDataset is None:
            print("[Annotate] WARNING: pybiomart not installed. Skipping BioMart annotation.", flush=True)
            print("  Install with: pip install pybiomart", flush=True)
            for col in ["SYMBOL", "pfam", "go_cc_id", "go_cc_term"]:
                df[col] = ""
            return df

        print("[Annotate] Connecting to Ensembl Plants BioMart ...", flush=True)
        try:
            dataset = BiomartDataset(
                name="athaliana_eg_gene",
                host="http://plants.ensembl.org",
                virtual_schema="plants_mart"
            )
        except Exception as e:
            print("[Annotate] WARNING: BioMart connection failed: {}".format(e), flush=True)
            for col in ["SYMBOL", "pfam", "go_cc_id", "go_cc_term"]:
                df[col] = ""
            return df

        gene_ids = df["jobs_new"].dropna().unique().tolist()
        print("[Annotate] Querying BioMart for {} gene IDs ...".format(len(gene_ids)), flush=True)

        try:
            annotations = dataset.query(
                attributes=[
                    "tair_locus", "external_gene_name",
                    "pfam", "go_id", "name_1006", "namespace_1003"
                ],
                filters={"link_ensembl_gene_id": gene_ids},
                use_attr_names=True
            )
        except Exception as e:
            print("[Annotate] WARNING: BioMart query failed: {}".format(e), flush=True)
            for col in ["SYMBOL", "pfam", "go_cc_id", "go_cc_term"]:
                df[col] = ""
            return df

        # GO Cellular Component
        go_cc = annotations[annotations["namespace_1003"] == "cellular_component"].copy()
        go_cc_summary = (
            go_cc.groupby("tair_locus", dropna=False)
            .agg(
                go_cc_id  =("go_id",    lambda x: ",".join(sorted(set(
                    v for v in x.dropna().astype(str) if v.strip())))),
                go_cc_term=("name_1006", lambda x: ",".join(sorted(set(
                    v for v in x.dropna().astype(str) if v.strip()))))
            )
            .reset_index()
        )

        def _first_non_na(s):
            vals = s.dropna().astype(str)
            vals = vals[vals.str.strip() != ""]
            return vals.iloc[0] if len(vals) > 0 else pd.NA

        annotation_summary = (
            annotations.groupby("tair_locus", dropna=False)
            .agg(
                SYMBOL=("external_gene_name", lambda x: ", ".join(sorted(set(
                    v for v in x.dropna().astype(str) if v.strip())))),
                pfam  =("pfam", _first_non_na)
            )
            .reset_index()
        )
        annotation_summary["SYMBOL"] = annotation_summary["SYMBOL"].fillna("")
        annotation_summary["pfam"]   = annotation_summary["pfam"].fillna("")

        df = df.merge(annotation_summary, left_on="jobs_new", right_on="tair_locus", how="left")
        df = df.drop(columns=["tair_locus"], errors="ignore")
        df = df.merge(go_cc_summary, left_on="jobs_new", right_on="tair_locus", how="left")
        df = df.drop(columns=["tair_locus"], errors="ignore")

        for col in ["SYMBOL", "pfam", "go_cc_id", "go_cc_term"]:
            df[col] = df[col].fillna("")

        print("[Annotate] BioMart annotation complete.", flush=True)
        return df

    def _annotate_interpro(self, df):
        df = df.copy()
        df["pfam"] = df["pfam"].replace("", pd.NA)

        pfam_summary = (
            df.assign(pfam=df["pfam"].fillna("Unknown"))
            .groupby("pfam", dropna=False)
            .agg(
                Frequency=("jobs_new", "size"),
                AT_ids   =("jobs_new", lambda x: ",".join(
                    pd.unique(x.dropna().astype(str))))
            )
            .reset_index()
            .sort_values("Frequency", ascending=False)
        )

        pfam_list = [p for p in pfam_summary["pfam"] if str(p).startswith("PF")]
        print("[Annotate] Querying InterPro for {} PFAM IDs ...".format(len(pfam_list)), flush=True)

        results = []
        for i, pfam_id in enumerate(pfam_list, 1):
            print("  [{}/{}] {}".format(i, len(pfam_list), pfam_id), flush=True)
            data = self._get_interpro_entry(pfam_id)
            if "error" not in data:
                meta      = data.get("metadata", {})
                name_dict = meta.get("name", {})
                results.append({
                    "pfam":       pfam_id,
                    "name":       name_dict.get("name", ""),
                    "short_name": name_dict.get("short", ""),
                    "type":       meta.get("type", "")
                })
            else:
                results.append({
                    "pfam":       pfam_id,
                    "name":       data["error"],
                    "short_name": "",
                    "type":       "unknown"
                })
            time.sleep(0.5)

        if results:
            results_df = pd.DataFrame(results)
        else:
            results_df = pd.DataFrame(columns=["pfam", "name", "short_name", "type"])

        pfam_final = pfam_summary.merge(results_df, on="pfam", how="left")
        pfam_final.to_csv(self.pfam_csv, index=False)
        print("[Annotate] PFAM summary saved: {}".format(self.pfam_csv), flush=True)

        if not results_df.empty:
            df = df.merge(results_df[["pfam", "name", "short_name", "type"]], on="pfam", how="left")
        else:
            df["name"] = df["short_name"] = df["type"] = ""

        for col in ["name", "short_name", "type"]:
            df[col] = df[col].fillna("")

        return df

    def _finalise_columns(self, df):
        cols_to_end = ["name", "short_name", "type", "go_cc_id", "go_cc_term"]
        existing    = [c for c in df.columns if c not in cols_to_end]
        present_end = [c for c in cols_to_end if c in df.columns]
        return df[existing + present_end]

    @staticmethod
    def _get_interpro_entry(pfam_id):
        url = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{}/".format(pfam_id)
        try:
            r = requests.get(url, headers={"Accept": "application/json"}, timeout=10)
            if r.status_code == 200:
                return r.json()
            return {"error": "Status {}".format(r.status_code)}
        except Exception as e:
            return {"error": str(e)}
