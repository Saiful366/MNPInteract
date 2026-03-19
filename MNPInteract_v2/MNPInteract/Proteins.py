"""Proteins.py — Protein list loading, UniProt sequence fetching, FASTA creation."""

import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests


def _wrap_fasta(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


class Proteins(object):

    def __init__(self, config):
        self.config = config
        self._session = requests.Session()
        self._sequence_cache = {}
        self.protein_rows = []

    def load(self):
        """Load protein list from the Stage 1 annotated CSV."""
        csv_path = self.config.csv_file
        if csv_path is None or not csv_path.exists():
            raise FileNotFoundError(
                "Annotated CSV not found: {}\n"
                "Run Stage 1 first: MNPInteract --S1".format(csv_path)
            )

        df = pd.read_csv(csv_path)
        for col in ("jobs", "jobs_new"):
            if col not in df.columns:
                raise ValueError(
                    "CSV '{}' must contain columns 'jobs' and 'jobs_new'.\n"
                    "Found: {}".format(csv_path, list(df.columns))
                )

        df = df[["jobs", "jobs_new"]].dropna().copy()
        df["jobs"]     = df["jobs"].astype(str).str.strip()
        df["jobs_new"] = df["jobs_new"].astype(str).str.strip()
        df = df.drop_duplicates()
        self.protein_rows = df.to_dict("records")
        print("[Proteins] Loaded {} protein entries.".format(len(self.protein_rows)), flush=True)

    def fetch_sequences(self, max_workers=8):
        unique_ids = sorted(set(r["jobs_new"] for r in self.protein_rows))
        print("\n[Proteins] Fetching sequences for {} gene IDs ...".format(len(unique_ids)), flush=True)

        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(self._fetch_one, gid): gid for gid in unique_ids}
            for fut in as_completed(futures):
                gid = futures[fut]
                try:
                    uid, seq = fut.result()
                    self._sequence_cache[gid] = (uid, seq)
                    print("  [OK] {} -> {}".format(gid, uid or "no UniProt"), flush=True)
                except Exception as e:
                    self._sequence_cache[gid] = (None, None)
                    print("  [FAIL] {}: {}".format(gid, e), flush=True)

    def _fetch_one(self, gene_id):
        uid = self._get_uniprot_id(gene_id)
        seq = self._get_sequence(uid) if uid else None
        return uid, seq

    def _get_uniprot_id(self, gene_id):
        if not gene_id.upper().startswith("AT"):
            return None
        url    = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query":  "(gene_exact:{}) AND (organism_id:3702)".format(gene_id),
            "format": "json",
            "fields": "accession",
            "size":   5,
        }
        r = self._session.get(url, params=params, timeout=20)
        r.raise_for_status()
        results = r.json().get("results", [])
        return results[0]["primaryAccession"] if results else None

    def _get_sequence(self, uniprot_id):
        url = "https://rest.uniprot.org/uniprotkb/{}.fasta".format(uniprot_id)
        r   = self._session.get(url, timeout=20)
        r.raise_for_status()
        lines = r.text.strip().splitlines()
        seq   = "".join(l.strip() for l in lines if not l.startswith(">"))
        return seq if seq else None

    def create_fasta_files(self):
        tasks_ready = []
        skipped     = []

        for row in self.protein_rows:
            jobs    = row["jobs"]
            gene_id = row["jobs_new"]
            folder  = self.config.pdb_dir / jobs

            if not folder.is_dir():
                skipped.append({"jobs": jobs, "jobs_new": gene_id, "reason": "directory not found"})
                print("[Skip] {} — folder not found".format(jobs), flush=True)
                continue

            uid, partner_seq = self._sequence_cache.get(gene_id, (None, None))
            if not partner_seq:
                skipped.append({"jobs": jobs, "jobs_new": gene_id, "reason": "sequence not found"})
                print("[Skip FASTA] No sequence for {}".format(gene_id), flush=True)
                continue

            fasta_path = folder / "input.fasta"
            with open(fasta_path, "w", newline="\n") as fh:
                fh.write(">PDLP5\n")
                fh.write(_wrap_fasta(self.config.pdlp5_seq) + "\n")
                fh.write(">{}\n".format(gene_id))
                fh.write(_wrap_fasta(partner_seq) + "\n")

            print("[FASTA] {}".format(fasta_path), flush=True)
            tasks_ready.append({
                "jobs":       jobs,
                "jobs_new":   gene_id,
                "target_dir": folder,
                "uniprot_id": uid or "",
            })

        print("\n[Proteins] FASTA ready: {}  |  Skipped: {}".format(
            len(tasks_ready), len(skipped)), flush=True)
        return tasks_ready, skipped
