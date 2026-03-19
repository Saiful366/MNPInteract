"""DeepTMHMM.py — Runs DeepTMHMM and parses topology results."""

import re
import csv
import shutil
import subprocess
import threading
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

MAX_RETRIES = 3
RETRY_DELAY = 30

_RE_LENGTH = re.compile(r"^#\s+(\S+)\s+Length:\s+\d+")
_RE_REGION = re.compile(r"^(\S+)\s+(signal|inside|outside|TMhelix|beta)\s+(\d+)\s+(\d+)")


def _join_ranges(ranges):
    return ",".join("{}-{}".format(s, e) for s, e in ranges)


class DeepTMHMM(object):

    def __init__(self, config):
        self.config = config
        self._report_lock = threading.Lock()
        self._run_report  = []
        self.run_report_csv  = config.output_dir / "deeptmhmm_run_report.csv"
        self.topo_report_csv = config.output_dir / "deeptmhmm_outside_inside_report.csv"

    def run(self, tasks):
        predict_py = self.config.deeptmhmm_dir / "predict.py"
        if not predict_py.exists():
            raise FileNotFoundError(
                "predict.py not found in Path_DeepTMHMM_dir:\n  {}\n"
                "Check that Path_DeepTMHMM_dir points to the unpacked "
                "DeepTMHMM-Academic-License directory.".format(predict_py)
            )

        pending = []
        for task in tasks:
            result_dir = task["target_dir"] / "result"
            if result_dir.exists() and any(result_dir.iterdir()):
                entry = self._make_entry(task, "skipped — result already exists")
                self._append_report(entry)
                print("[Skip DeepTMHMM] {} — already done".format(task["jobs"]), flush=True)
            else:
                pending.append(task)

        print("\n[DeepTMHMM] Jobs to run: {}  |  Workers: {}\n".format(
            len(pending), self.config.max_workers_gpu), flush=True)

        with ThreadPoolExecutor(max_workers=self.config.max_workers_gpu) as executor:
            futures = {executor.submit(self._run_one, task, predict_py): task
                       for task in pending}
            for i, fut in enumerate(as_completed(futures), 1):
                task = futures[fut]
                try:
                    entry = fut.result()
                except Exception as e:
                    entry = self._make_entry(task, "unexpected error: {}".format(e))
                self._append_report(entry)
                print("[DeepTMHMM {}/{}] {} — {}".format(
                    i, len(pending), task["jobs"], entry["status"]), flush=True)

        return self._run_report

    def parse_results(self):
        print("\n[DeepTMHMM] Parsing topology results ...", flush=True)
        pdb_dir = self.config.pdb_dir
        rows    = []
        skipped = []

        for subdir in sorted(pdb_dir.iterdir()):
            if not subdir.is_dir():
                continue
            md_file = subdir / "result" / "deeptmhmm_results.md"
            if not md_file.exists():
                print("  [Missing] {} — no deeptmhmm_results.md".format(subdir.name), flush=True)
                skipped.append(subdir.name)
                continue

            proteins_order, regions = self._parse_md(md_file)
            if not proteins_order:
                print("  [Empty] {} — no proteins parsed".format(subdir.name), flush=True)
                skipped.append(subdir.name)
                continue

            pdlp5   = proteins_order[0]
            partner = proteins_order[1] if len(proteins_order) > 1 else ""

            pdlp5_outside  = _join_ranges(regions[pdlp5]["outside"])
            partner_inside = _join_ranges(regions.get(partner, {}).get("inside", [])) or "NA"

            rows.append({
                "folder":                subdir.name,
                "partner_protein":       partner,
                "PDLP5_outside_range":   pdlp5_outside,
                "partner_inside_ranges": partner_inside,
            })

        df = pd.DataFrame(rows)
        df.to_csv(self.topo_report_csv, index=False, quoting=csv.QUOTE_ALL)
        print("[DeepTMHMM] Topology CSV: {}".format(self.topo_report_csv), flush=True)
        print("  Parsed: {}  |  Skipped: {}".format(len(rows), len(skipped)), flush=True)
        return self.topo_report_csv

    def _run_one(self, task, predict_py):
        target_dir = task["target_dir"]
        result_dir = target_dir / "result"
        fasta_path = target_dir / "input.fasta"

        if result_dir.exists():
            shutil.rmtree(result_dir)

        for attempt in range(1, MAX_RETRIES + 1):
            try:
                cmd = [
                    str(self.config.deeptmhmm_py),
                    str(predict_py),
                    "--fasta",      str(fasta_path),
                    "--output-dir", str(result_dir),
                ]
                proc = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=600,
                    cwd=str(self.config.deeptmhmm_dir),
                )
                if proc.returncode != 0:
                    raise RuntimeError(
                        "predict.py exited {}\n{}".format(
                            proc.returncode, proc.stderr[-500:])
                    )
                return self._make_entry(task, "completed")
            except Exception as e:
                if attempt < MAX_RETRIES:
                    print("  [Retry {}/{}] {}: {}".format(
                        attempt, MAX_RETRIES, task["jobs"], e), flush=True)
                    time.sleep(RETRY_DELAY)
                else:
                    return self._make_entry(
                        task, "failed after {} attempts: {}".format(MAX_RETRIES, e))

    @staticmethod
    def _parse_md(md_file):
        try:
            with open(md_file, "r", encoding="utf-8-sig", errors="replace") as fh:
                lines = [l.strip() for l in fh if l.strip()]
        except Exception as e:
            print("  [Error reading] {}: {}".format(md_file, e), flush=True)
            return [], {}

        proteins_order = []
        regions        = {}

        for line in lines:
            m = _RE_LENGTH.match(line)
            if m:
                prot = m.group(1)
                if prot not in proteins_order:
                    proteins_order.append(prot)
                regions.setdefault(prot, {"inside": [], "outside": []})
                continue
            m = _RE_REGION.match(line)
            if m:
                prot, rtype, start, end = m.groups()
                start, end = int(start), int(end)
                regions.setdefault(prot, {"inside": [], "outside": []})
                if rtype == "inside":
                    regions[prot]["inside"].append((start, end))
                elif rtype == "outside":
                    regions[prot]["outside"].append((start, end))

        return proteins_order, regions

    def _append_report(self, entry):
        with self._report_lock:
            self._run_report.append(entry)
            pd.DataFrame(self._run_report).to_csv(self.run_report_csv, index=False)

    @staticmethod
    def _make_entry(task, status):
        return {
            "jobs":              task["jobs"],
            "jobs_new":          task["jobs_new"],
            "matched_directory": str(task["target_dir"]),
            "uniprot_id":        task.get("uniprot_id", ""),
            "status":            status,
        }
