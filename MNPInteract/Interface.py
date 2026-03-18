"""Interface.py — Atom-level interface detection with PAE filtering."""

import gzip
import json
import pickle
import sys
import threading
import types
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser


# ── JAX compatibility shim ────────────────────────────────────────────────────

def _extract_array_like(obj):
    if isinstance(obj, np.ndarray):
        return obj
    if isinstance(obj, (list, tuple)):
        for item in obj:
            arr = _extract_array_like(item)
            if arr is not None:
                return arr
    if hasattr(obj, "shape") and hasattr(obj, "dtype"):
        try:
            return np.asarray(obj)
        except Exception:
            pass
    return None


def reconstruct_device_array(*args, **kwargs):
    for obj in args:
        arr = _extract_array_like(obj)
        if arr is not None:
            return np.asarray(arr)
    return np.array([])


_fake_jax = types.ModuleType("jax._src.device_array")


class DeviceArray(np.ndarray):
    pass


_fake_jax.DeviceArray = DeviceArray
_fake_jax.reconstruct_device_array = reconstruct_device_array
sys.modules["jax._src.device_array"] = _fake_jax


class CompatUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "jax._src.device_array":
            if name == "DeviceArray":
                return DeviceArray
            if name == "reconstruct_device_array":
                return reconstruct_device_array
        return super(CompatUnpickler, self).find_class(module, name)


# ── Module-level constants (set by worker initialiser) ────────────────────────

_CHAIN1          = "B"
_CHAIN2          = "C"
_CONTACT_ATOMS   = {"C", "CA", "CB"}
_DISTANCE_CUTOFF = 6.0
_PAE_CUTOFF      = 25.0


def _worker_init(chain1, chain2, contact_atoms, dist_cutoff, pae_cutoff):
    global _CHAIN1, _CHAIN2, _CONTACT_ATOMS, _DISTANCE_CUTOFF, _PAE_CUTOFF
    _CHAIN1          = chain1
    _CHAIN2          = chain2
    _CONTACT_ATOMS   = contact_atoms
    _DISTANCE_CUTOFF = dist_cutoff
    _PAE_CUTOFF      = pae_cutoff


# ── Per-folder helpers ────────────────────────────────────────────────────────

def _load_pickle(pkl_file):
    if str(pkl_file).endswith(".gz"):
        with open(pkl_file, "rb") as fh:
            return CompatUnpickler(gzip.open(fh)).load()
    with open(pkl_file, "rb") as fh:
        return CompatUnpickler(fh).load()


def _is_standard(residue):
    return residue.id[0] == " "


def _residue_label(chain_id, residue):
    return "{}:{} {}".format(chain_id, residue.get_resname(), residue.id[1])


def _build_res_index_map(model):
    mapping = {}
    idx = 0
    for chain in model:
        for res in chain:
            if _is_standard(res):
                mapping[(chain.id, res.id)] = idx
                idx += 1
    return mapping


def _min_atom_dist(res1, res2):
    best = None
    for a1 in res1:
        if a1.get_id() not in _CONTACT_ATOMS:
            continue
        for a2 in res2:
            if a2.get_id() not in _CONTACT_ATOMS:
                continue
            d = float(a1 - a2)
            if best is None or d < best:
                best = d
    return best


def _find_ranked_pdb(subdir):
    pdb = subdir / "ranked_0.pdb"
    if pdb.exists():
        return pdb
    pdbs = sorted(subdir.glob("*.pdb"))
    return pdbs[0] if pdbs else None


def _find_result_pkl(subdir):
    ranking_json = subdir / "ranking_debug.json"
    if ranking_json.exists():
        try:
            with open(ranking_json) as fh:
                ranking = json.load(fh)
            top_model = ranking["order"][0]
            for ext in (".pkl", ".pkl.gz"):
                candidate = subdir / "result_{}{}".format(top_model, ext)
                if candidate.exists():
                    return candidate
            print("  [WARNING] PKL for top model '{}' not found in {}. Falling back.".format(
                top_model, subdir.name), flush=True)
        except Exception as e:
            print("  [WARNING] ranking_debug.json read error: {}".format(e), flush=True)

    for pattern in ("result_model*.pkl.gz", "result_model*.pkl", "*.pkl.gz", "*.pkl"):
        matches = sorted(subdir.glob(pattern))
        if matches:
            return matches[0]
    return None


def _save_outputs(subdir, df, df_unique, iface1, iface2):
    df.to_csv(subdir / "BC_res_int_atom.csv", index=False)
    df_unique.to_csv(subdir / "BC_res_int_atom_unique_pairs.csv", index=False)
    (subdir / "BC_res_int_atom_{}_unique_residues.txt".format(_CHAIN1)).write_text(
        ",".join(map(str, iface1)) + "\n", encoding="utf-8")
    (subdir / "BC_res_int_atom_{}_unique_residues.txt".format(_CHAIN2)).write_text(
        ",".join(map(str, iface2)) + "\n", encoding="utf-8")
    pymol1 = "select interface_{}, chain {} and resi {}".format(
        _CHAIN1, _CHAIN1, "+".join(map(str, iface1)))
    pymol2 = "select interface_{}, chain {} and resi {}".format(
        _CHAIN2, _CHAIN2, "+".join(map(str, iface2)))
    (subdir / "BC_res_int_atom_pymol_selections.txt").write_text(
        pymol1 + "\n" + pymol2 + "\n", encoding="utf-8")


def _save_empty_outputs(subdir):
    for fname in (
        "BC_res_int_atom.csv",
        "BC_res_int_atom_unique_pairs.csv",
        "BC_res_int_atom_{}_unique_residues.txt".format(_CHAIN1),
        "BC_res_int_atom_{}_unique_residues.txt".format(_CHAIN2),
        "BC_res_int_atom_pymol_selections.txt",
    ):
        (subdir / fname).write_text("", encoding="utf-8")


def _report_entry(subdir, pdb, pkl, status, n_pairs):
    return {
        "folder":          subdir.name,
        "pdb_file":        str(pdb),
        "pkl_file":        str(pkl),
        "status":          status,
        "interface_pairs": n_pairs,
        "pkl_error":       "pickle load error" in str(status),
    }


def process_folder(subdir):
    """Entry point for each subprocess worker."""
    pdb_file = _find_ranked_pdb(subdir)
    pkl_file = _find_result_pkl(subdir)

    if pdb_file is None:
        return _report_entry(subdir, "", "", "ranked_0.pdb not found", 0)
    if pkl_file is None:
        return _report_entry(subdir, pdb_file, "", ".pkl file not found", 0)

    print("\n" + "-"*60, flush=True)
    print("Folder : {}".format(subdir.name), flush=True)
    print("PDB    : {}".format(pdb_file.name), flush=True)
    print("PKL    : {}".format(pkl_file.name), flush=True)

    try:
        try:
            pickle_dict = _load_pickle(pkl_file)
        except Exception as e:
            print("  ERROR loading pickle: {}".format(e), flush=True)
            return _report_entry(subdir, pdb_file, pkl_file,
                                 "pickle load error: {}".format(e), 0)

        print("Pickle keys: {}".format(list(pickle_dict.keys())), flush=True)

        if "predicted_aligned_error" not in pickle_dict:
            return _report_entry(subdir, pdb_file, pkl_file,
                                 "predicted_aligned_error not found in pickle", 0)

        pae = np.array(pickle_dict["predicted_aligned_error"], dtype=float)
        print("PAE shape: {}".format(pae.shape), flush=True)

        parser = PDBParser(QUIET=True)
        with open(pdb_file, "r", encoding="utf-8", errors="replace") as fh:
            structure = parser.get_structure("model", fh)
        model = structure[0]

        for chain_id in (_CHAIN1, _CHAIN2):
            if chain_id not in model:
                return _report_entry(subdir, pdb_file, pkl_file,
                                     "Chain {} not found in PDB".format(chain_id), 0)

        res_index = _build_res_index_map(model)
        n_res     = len(res_index)

        if pae.shape != (n_res, n_res):
            return _report_entry(
                subdir, pdb_file, pkl_file,
                "PAE shape {} != residue count ({},{}). "
                "PDB and PKL may be from different runs.".format(
                    pae.shape, n_res, n_res), 0)

        chain1_res = [r for r in model[_CHAIN1] if _is_standard(r)]
        chain2_res = [r for r in model[_CHAIN2] if _is_standard(r)]

        print("Chain {}: {} residues | Chain {}: {} residues".format(
            _CHAIN1, len(chain1_res), _CHAIN2, len(chain2_res)), flush=True)
        print("Cutoffs: distance <= {} A, PAE < {}".format(
            _DISTANCE_CUTOFF, _PAE_CUTOFF), flush=True)

        rows = []
        for r1 in chain1_res:
            idx1 = res_index[(_CHAIN1, r1.id)]
            for r2 in chain2_res:
                idx2    = res_index[(_CHAIN2, r2.id)]
                dist    = _min_atom_dist(r1, r2)
                if dist is None or dist > _DISTANCE_CUTOFF:
                    continue
                pae_val = float(pae[idx2][idx1])
                if pae_val >= _PAE_CUTOFF:
                    continue
                rows.append({
                    _CHAIN1:      _residue_label(_CHAIN1, r1),
                    _CHAIN2:      _residue_label(_CHAIN2, r2),
                    "Distance_A": round(dist, 7),
                    "PAE":        round(pae_val, 7),
                })

        if not rows:
            print("No interface pairs found with current cutoffs.", flush=True)
            _save_empty_outputs(subdir)
            return _report_entry(subdir, pdb_file, pkl_file, "no interface pairs found", 0)

        df = (
            pd.DataFrame(rows)
            .sort_values(["Distance_A", "PAE"])
            .reset_index(drop=True)
        )
        print("Found {} interface pairs".format(len(df)), flush=True)

        iface1 = sorted(set(int(x.split()[-1]) for x in df[_CHAIN1]))
        iface2 = sorted(set(int(x.split()[-1]) for x in df[_CHAIN2]))

        seen, unique_rows = set(), []
        for _, row in df.iterrows():
            key = frozenset([row[_CHAIN1], row[_CHAIN2]])
            if key not in seen:
                seen.add(key)
                unique_rows.append(row)
        df_unique = pd.DataFrame(unique_rows).reset_index(drop=True)

        _save_outputs(subdir, df, df_unique, iface1, iface2)
        return _report_entry(subdir, pdb_file, pkl_file, "success", len(df))

    except Exception as exc:
        import traceback
        traceback.print_exc()
        return _report_entry(subdir, pdb_file, pkl_file, "error: {}".format(exc), 0)


# ── Interface class ───────────────────────────────────────────────────────────

class Interface(object):

    CHAIN1 = "B"
    CHAIN2 = "C"

    def __init__(self, config):
        self.config     = config
        self.report_csv = config.output_dir / "interface_processing_report_atom.csv"
        self._summary   = []
        self._lock      = threading.Lock()

    def run(self):
        pdb_dir = self.config.pdb_dir
        subdirs = [d for d in sorted(pdb_dir.iterdir()) if d.is_dir()]

        if not subdirs:
            print("[Interface] No subdirectories found in pdb_dir.", flush=True)
            return []

        print("\n" + "="*60, flush=True)
        print("Stage 3 — Interface detection (atom-level)", flush=True)
        print("="*60, flush=True)
        print("Folders : {}".format(len(subdirs)), flush=True)
        print("Cutoffs : distance <= {} A, PAE < {}".format(
            self.config.distance_cutoff, self.config.pae_cutoff), flush=True)
        print("Workers : {}".format(self.config.max_workers_cpu), flush=True)

        pending = []
        for sd in subdirs:
            pairs_csv = sd / "BC_res_int_atom_unique_pairs.csv"
            if pairs_csv.exists() and pairs_csv.stat().st_size > 0:
                print("[Skip] {} — already processed".format(sd.name), flush=True)
                self._summary.append({
                    "folder": sd.name, "pdb_file": "", "pkl_file": "",
                    "status": "skipped — already processed",
                    "interface_pairs": "", "pkl_error": False,
                })
            else:
                pending.append(sd)

        print("\nPending: {} folders\n".format(len(pending)), flush=True)

        init_args = (
            self.CHAIN1, self.CHAIN2,
            _CONTACT_ATOMS,
            self.config.distance_cutoff,
            self.config.pae_cutoff,
        )

        with ProcessPoolExecutor(
            max_workers=self.config.max_workers_cpu,
            initializer=_worker_init,
            initargs=init_args,
        ) as executor:
            futures = {executor.submit(process_folder, sd): sd for sd in pending}
            for i, fut in enumerate(as_completed(futures), 1):
                sd = futures[fut]
                try:
                    result = fut.result()
                except Exception as e:
                    result = _report_entry(sd, "", "", "unexpected error: {}".format(e), 0)
                with self._lock:
                    self._summary.append(result)
                print("[Interface {}/{}] {}: {} ({} pairs)".format(
                    i, len(pending), sd.name,
                    result["status"], result["interface_pairs"]), flush=True)

        df = pd.DataFrame(self._summary)
        df.to_csv(self.report_csv, index=False)
        print("\n[Interface] Report saved: {}".format(self.report_csv), flush=True)
        return self._summary
