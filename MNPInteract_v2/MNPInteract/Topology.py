"""Topology.py — Topology-based LIKELY/UNLIKELY verdict assignment."""

from pathlib import Path
import pandas as pd


def _parse_ranges(range_str):
    s = str(range_str).strip()
    if not s or s.lower() in ("nan", "na", "n/a", ""):
        return []
    result = []
    for seg in s.replace(";", ",").split(","):
        seg = seg.strip()
        if "-" in seg:
            parts = seg.split("-")
            try:
                result.append((int(parts[0]), int(parts[1])))
            except ValueError:
                pass
    return result


def _in_any_range(res_num, ranges):
    return any(start <= res_num <= end for start, end in ranges)


def _extract_res_num(label):
    return int(str(label).strip().split()[-1])


class Topology(object):

    def __init__(self, config, topo_csv):
        self.config     = config
        self.topo_csv   = topo_csv
        self.report_csv = config.output_dir / "topology_contact_report.csv"

    def run(self):
        print("\n" + "="*60, flush=True)
        print("Stage 4 — Topology filter", flush=True)
        print("="*60, flush=True)

        if not self.topo_csv.exists():
            print("[Topology] WARNING: topology file not found: {}".format(self.topo_csv), flush=True)
            print("           Skipping topology filter.", flush=True)
            return pd.DataFrame()

        topo_df = pd.read_csv(self.topo_csv)
        print("Loaded: {}  ({} rows)".format(self.topo_csv.name, len(topo_df)), flush=True)

        required = {"folder", "partner_protein", "PDLP5_outside_range", "partner_inside_ranges"}
        missing  = required - set(topo_df.columns)
        if missing:
            raise ValueError("Topology CSV missing columns: {}".format(missing))

        rows = []
        for _, trow in topo_df.iterrows():
            row = self._process_row(trow)
            rows.append(row)
            print("  {:<50s}  verdict: {}".format(
                row["folder"], row["interaction_verdict"]), flush=True)

        report_df = pd.DataFrame(rows)
        report_df.to_csv(self.report_csv, index=False)
        print("\n[Topology] Report saved: {}".format(self.report_csv), flush=True)

        counts = report_df["interaction_verdict"].value_counts()
        print("\nVerdict summary:\n{}".format(counts.to_string()), flush=True)

        return report_df

    def _process_row(self, trow):
        folder_name    = str(trow["folder"]).strip()
        partner        = str(trow["partner_protein"]).strip()
        pdlp5_outside  = _parse_ranges(trow["PDLP5_outside_range"])
        partner_inside = _parse_ranges(trow["partner_inside_ranges"])

        base = {
            "folder":                folder_name,
            "partner_protein":       partner,
            "PDLP5_outside_range":   trow["PDLP5_outside_range"],
            "partner_inside_ranges": trow["partner_inside_ranges"],
        }

        if not partner_inside:
            return dict(list(base.items()) + [
                ("total_pairs", 0), ("unlikely_pairs_count", 0),
                ("unlikely_pairs", ""), ("interaction_verdict", "LIKELY")])

        folder_path = self.config.pdb_dir / folder_name
        if not folder_path.is_dir():
            return dict(list(base.items()) + [
                ("total_pairs", 0), ("unlikely_pairs_count", 0),
                ("unlikely_pairs", ""), ("interaction_verdict", "FOLDER_NOT_FOUND")])

        pairs_csv = folder_path / "BC_res_int_atom_unique_pairs.csv"
        if not pairs_csv.exists():
            return dict(list(base.items()) + [
                ("total_pairs", 0), ("unlikely_pairs_count", 0),
                ("unlikely_pairs", ""), ("interaction_verdict", "MISSING_CSV")])

        try:
            pairs_df = pd.read_csv(pairs_csv)
        except pd.errors.EmptyDataError:
            return dict(list(base.items()) + [
                ("total_pairs", 0), ("unlikely_pairs_count", 0),
                ("unlikely_pairs", ""), ("interaction_verdict", "NO_PAIRS")])

        if pairs_df.empty:
            return dict(list(base.items()) + [
                ("total_pairs", 0), ("unlikely_pairs_count", 0),
                ("unlikely_pairs", ""), ("interaction_verdict", "NO_PAIRS")])

        b_cols = [c for c in pairs_df.columns if c.strip().upper() == "B"]
        c_cols = [c for c in pairs_df.columns if c.strip().upper() == "C"]
        if not b_cols or not c_cols:
            return dict(list(base.items()) + [
                ("total_pairs", len(pairs_df)), ("unlikely_pairs_count", 0),
                ("unlikely_pairs", ""), ("interaction_verdict", "COLUMN_ERROR")])

        b_col = b_cols[0]
        c_col = c_cols[0]
        total_pairs    = len(pairs_df)
        unlikely_pairs = []

        for _, pair_row in pairs_df.iterrows():
            try:
                b_num = _extract_res_num(str(pair_row[b_col]))
                c_num = _extract_res_num(str(pair_row[c_col]))
            except Exception:
                continue
            if _in_any_range(b_num, pdlp5_outside) and _in_any_range(c_num, partner_inside):
                unlikely_pairs.append("{} <-> {}".format(
                    str(pair_row[b_col]).strip(), str(pair_row[c_col]).strip()))

        n_unlikely = len(unlikely_pairs)
        verdict    = "UNLIKELY" if n_unlikely > 0 else "LIKELY"

        return dict(list(base.items()) + [
            ("total_pairs", total_pairs),
            ("unlikely_pairs_count", n_unlikely),
            ("unlikely_pairs", " ; ".join(unlikely_pairs)),
            ("interaction_verdict", verdict),
        ])
