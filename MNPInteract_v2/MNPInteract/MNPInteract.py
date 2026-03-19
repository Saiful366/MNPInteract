"""MNPInteract.py — Main entry point and pipeline orchestrator.

All steps run on Nova (HPC).

Step 1  Login node  (internet required, no GPU):
    MNPInteract --S1

Step 2  SLURM GPU job:
    MNPInteract --S2-4

Step 3  Login node  (no GPU needed):
    MNPInteract --S5
"""

import argparse
import sys
import logging
from pathlib import Path

log_filename = "./MNPInteract.log"
logging.basicConfig(
    filename=log_filename, filemode="w", level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


class _TeeLogger(object):
    def __init__(self, log_file):
        self._terminal = sys.stdout
        self._log      = open(log_file, "a")

    def write(self, msg):
        self._terminal.write(msg)
        self._log.write(msg)

    def flush(self):
        self._terminal.flush()
        self._log.flush()


sys.stdout = _TeeLogger(log_filename)
sys.stderr = _TeeLogger(log_filename)

from .Config     import Config
from .Annotate   import Annotate
from .Proteins   import Proteins
from .DeepTMHMM  import DeepTMHMM
from .Interface  import Interface
from .Topology   import Topology
from .Confidence import Confidence


CONF_TEMPLATE = """\
# MNPInteract configuration file
# Lines starting with '#' are comments.
# Format:  Key : Value
# All paths must be absolute paths on Nova (Lustre filesystem).

# ── INPUT ──────────────────────────────────────────────────────────────────────

# REQUIRED for --S1 — AlphaPulldown scoring output
Path_interpae_csv     : /lustre/hdd/LAS/<your-lab>/<username>/alphapulldown/predictions_with_good_interpae.csv

# REQUIRED for --S2-4 — directory containing AlphaPulldown prediction sub-folders
Path_PDB_dir          : /lustre/hdd/LAS/<your-lab>/<username>/alphapulldown

# OPTIONAL — PDLP5 FASTA. Leave blank to use the built-in AT2G33330 sequence.
Path_PDLP5_fasta      :

# REQUIRED — all output files go here (created automatically if absent)
Path_output           : /lustre/hdd/LAS/<your-lab>/<username>/alphapulldown/MNPInteract_output

# ── DEEPTMHMM (--S2-4 only) ────────────────────────────────────────────────────

# Path to the unpacked DeepTMHMM directory (must directly contain predict.py)
Path_DeepTMHMM_dir    : /lustre/hdd/LAS/<your-lab>/<username>/alphapulldown/DeepTMHMM-Academic-License-v1.0

# Python binary inside the DeepTMHMM virtual environment
Path_DeepTMHMM_venv   : /lustre/hdd/LAS/<your-lab>/<username>/alphapulldown/deeptmhmm-venv/bin/python

# ── PARALLELISM ────────────────────────────────────────────────────────────────

Max_workers_gpu       : 4
Max_workers_cpu       : 16

# ── CONTACT CUTOFFS ────────────────────────────────────────────────────────────

Distance_cutoff       : 6.0
PAE_cutoff            : 25.0

# ── SCORE CUTOFFS (--S1 filtering, OR logic) ───────────────────────────────────

IPTM_cutoff           : 0.30
DockQ_cutoff          : 0.23
Piscore_cutoff        : 0
"""


def _build_parser():
    p = argparse.ArgumentParser(
        prog="MNPInteract",
        description=(
            "Post-AlphaPulldown pipeline for identifying high-confidence\n"
            "PDLP5-interacting proteins. All steps run on Nova (HPC).\n\n"
            "Step 1 - login node (internet required, no GPU):\n"
            "    MNPInteract --S1\n\n"
            "Step 2 - SLURM GPU job:\n"
            "    MNPInteract --S2-4\n\n"
            "Step 3 - login node (no GPU):\n"
            "    MNPInteract --S5"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--S1", action="store_true",
        help="Stage 1: score filtering + PFAM/GO annotation  (login node)"
    )
    group.add_argument(
        "--S2-4", dest="S2_4", action="store_true",
        help="Stages 2-4: DeepTMHMM + interface + topology   (SLURM GPU job)"
    )
    group.add_argument(
        "--S5", action="store_true",
        help="Stage 5: final high-confidence interactors      (login node)"
    )
    group.add_argument(
        "--print-template", action="store_true",
        help="Print a conf.txt template and exit"
    )

    p.add_argument(
        "--conf", default="conf.txt",
        help="Path to conf.txt (default: ./conf.txt)"
    )
    return p


def _run_stage1(config):
    print("\n" + "="*60, flush=True)
    print("STAGE 1 — Score filtering and PFAM/GO annotation", flush=True)
    print("="*60, flush=True)

    annotator     = Annotate(config)
    annotated_csv = annotator.run()

    print("\n" + "-"*60, flush=True)
    print("  Stage 1 complete.", flush=True)
    print("  Output: {}".format(annotated_csv), flush=True)
    print("-"*60, flush=True)
    print("  NEXT STEP:", flush=True)
    print("  Submit a SLURM GPU job and run:", flush=True)
    print("    MNPInteract --S2-4", flush=True)
    print("-"*60 + "\n", flush=True)


def _run_stage24(config):
    print("\n" + "="*60, flush=True)
    print("STAGES 2-4 — DeepTMHMM -> Interface -> Topology", flush=True)
    print("="*60, flush=True)

    annotated_csv = (
        config.output_dir /
        "Prioritized candidate list_with_PFAM_GO_annotations.csv"
    )
    if not annotated_csv.exists():
        raise FileNotFoundError(
            "Stage 1 output not found:\n  {}\n"
            "Run Stage 1 first on the login node:\n"
            "  MNPInteract --S1".format(annotated_csv)
        )

    config.csv_file = annotated_csv

    print("\n[Prep] Loading protein list from Stage 1 output ...", flush=True)
    proteins = Proteins(config)
    proteins.load()
    proteins.fetch_sequences(max_workers=8)
    tasks_ready, skipped = proteins.create_fasta_files()
    if skipped:
        print("  [{} entries skipped]".format(len(skipped)), flush=True)

    # Stage 2
    print("\n" + "-"*60, flush=True)
    print("Stage 2 — DeepTMHMM topology prediction", flush=True)
    dtm          = DeepTMHMM(config)
    dtm.run(tasks_ready)
    dtm_topo_csv = dtm.parse_results()

    # Stage 3 — only process folders from the filtered candidate list
    print("\n" + "-"*60, flush=True)
    print("Stage 3 — Atom-level interface detection", flush=True)
    folder_names = [row["jobs"] for row in proteins.protein_rows]
    iface = Interface(config)
    iface.run(folder_names=folder_names)

    # Stage 4
    print("\n" + "-"*60, flush=True)
    print("Stage 4 — Topology filter", flush=True)
    if dtm_topo_csv.exists():
        topo = Topology(config, dtm_topo_csv)
        topo.run()
    else:
        print("  WARNING: DeepTMHMM topology CSV not found. Stage 4 skipped.", flush=True)

    topo_csv = config.output_dir / "topology_contact_report.csv"
    print("\n" + "-"*60, flush=True)
    print("  Stages 2-4 complete.", flush=True)
    print("  Output: {}".format(topo_csv), flush=True)
    print("-"*60, flush=True)
    print("  NEXT STEP:", flush=True)
    print("  On the login node run:", flush=True)
    print("    MNPInteract --S5", flush=True)
    print("-"*60 + "\n", flush=True)


def _run_stage5(config):
    print("\n" + "="*60, flush=True)
    print("STAGE 5 — Final high-confidence interactor generation", flush=True)
    print("="*60, flush=True)

    confidence = Confidence(config)
    final_df   = confidence.run()

    final_path = config.output_dir / "Final_High_Confidence_Interactors_only.csv"
    print("\n" + "-"*60, flush=True)
    print("  Stage 5 complete.", flush=True)
    print("  Final output : {}".format(final_path), flush=True)
    print("  Interactors  : {}".format(len(final_df)), flush=True)
    print("-"*60 + "\n", flush=True)


def main():
    parser = _build_parser()
    args   = parser.parse_args()

    if args.print_template:
        print(CONF_TEMPLATE)
        sys.exit(0)

    if args.S1:
        stage       = "1"
        stage_label = "--S1"
    elif args.S2_4:
        stage       = "2-4"
        stage_label = "--S2-4"
    else:
        stage       = "5"
        stage_label = "--S5"

    print("="*60, flush=True)
    print("  MNPInteract — PDLP5 interaction confidence pipeline", flush=True)
    print("  Running: {}".format(stage_label), flush=True)
    print("="*60, flush=True)

    print("\n[Config] Loading configuration ...", flush=True)
    config = Config(conf_file=args.conf, stage=stage)
    print(config, flush=True)

    if stage == "1":
        _run_stage1(config)
    elif stage == "2-4":
        _run_stage24(config)
    elif stage == "5":
        _run_stage5(config)

    print("="*60, flush=True)
    print("  MNPInteract — Done.", flush=True)
    print("  Output directory : {}".format(config.output_dir), flush=True)
    print("  Log file         : {}".format(Path(log_filename).resolve()), flush=True)
    print("="*60, flush=True)


if __name__ == "__main__":
    main()
