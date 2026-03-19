"""Config.py — Stage-aware configuration reader for MNPInteract.

Stage 1   needs: Path_interpae_csv, Path_output
Stage 2-4 needs: Path_PDB_dir, Path_DeepTMHMM_dir, Path_DeepTMHMM_venv, Path_output
Stage 5   needs: Path_output only
"""

from pathlib import Path

PDLP5_BUILTIN_SEQ = (
    "MIKTKTTSLLCFLLTAVILMNPSSSSPTDNYIYAVCSPAKFSPSSGYETNLNSLLSSFVT"
    "STAQTRYANFTVPTGKPEPTVTVYGIYQCRGDLDPTACSTCVSSAVAQVGALCSNSYSGF"
    "LQMENCLIRYDNKSFLGVQDKTLILNKCGQPMEFNDQDALTKASDVIGSLGTGDGSYRTG"
    "GNGNVQGVAQCSGDLSTSQCQDCLSDAIGRLKSDCGMAQGGYVYLSKCYARFSVGGSHAR"
    "QTPGPNFGHEGEKGNKDDNGVGKTLAIIIGIVTLIILLVVFLAFVGKCCRKLQDEKWCK"
)

DEFAULTS = {
    "Max_workers_gpu": 4,
    "Max_workers_cpu": 16,
    "Distance_cutoff": 6.0,
    "PAE_cutoff":      25.0,
    "IPTM_cutoff":     0.30,
    "DockQ_cutoff":    0.23,
    "Piscore_cutoff":  0.0,
}


class Config(object):

    def __init__(self, conf_file="conf.txt", stage="1"):
        self._conf_file = Path(conf_file)
        self._stage     = stage
        self._raw       = {}
        self._parse_conf()
        self._resolve()

    def _parse_conf(self):
        if not self._conf_file.exists():
            raise FileNotFoundError(
                "conf.txt not found: {}\n"
                "Run 'MNPInteract --print-template' to generate a template.".format(
                    self._conf_file)
            )
        with open(self._conf_file) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or ":" not in line:
                    continue
                key, _, value = line.partition(":")
                self._raw[key.strip()] = value.strip()

    def _resolve(self):
        def get(key, required=True, default=None):
            val = self._raw.get(key, "").strip()
            if not val:
                if required:
                    raise ValueError(
                        "conf.txt: required key '{}' is missing or empty.\n"
                        "Run 'MNPInteract --print-template' to see all keys.".format(key)
                    )
                return default
            return val

        s = self._stage

        # Always required
        self.output_dir = Path(get("Path_output"))

        # Stage 1
        need_1 = (s == "1")
        raw_ip = get("Path_interpae_csv", required=need_1, default="")
        self.interpae_csv = Path(raw_ip) if raw_ip else None

        # Stage 2-4
        need_24 = (s == "2-4")
        raw_pdb  = get("Path_PDB_dir",        required=need_24, default="")
        raw_dtm  = get("Path_DeepTMHMM_dir",  required=need_24, default="")
        raw_venv = get("Path_DeepTMHMM_venv", required=need_24, default="")

        self.pdb_dir       = Path(raw_pdb)  if raw_pdb  else None
        self.deeptmhmm_dir = Path(raw_dtm)  if raw_dtm  else None
        self.deeptmhmm_py  = Path(raw_venv) if raw_venv else None

        # Runtime attributes set by orchestrator
        self.csv_file  = None
        self.test_file = None

        # PDLP5 sequence
        raw_fasta = get("Path_PDLP5_fasta", required=False, default="")
        self.pdlp5_fasta = Path(raw_fasta) if raw_fasta else None
        if self.pdlp5_fasta and self.pdlp5_fasta.exists():
            self.pdlp5_seq = self._read_fasta(self.pdlp5_fasta)
        else:
            self.pdlp5_seq = PDLP5_BUILTIN_SEQ

        # Numeric parameters
        def fget(k):
            return float(get(k, required=False, default=str(DEFAULTS[k])))
        def iget(k):
            return int(get(k, required=False, default=str(DEFAULTS[k])))

        self.max_workers_gpu = iget("Max_workers_gpu")
        self.max_workers_cpu = iget("Max_workers_cpu")
        self.distance_cutoff = fget("Distance_cutoff")
        self.pae_cutoff      = fget("PAE_cutoff")
        self.iptm_cutoff     = fget("IPTM_cutoff")
        self.dockq_cutoff    = fget("DockQ_cutoff")
        self.piscore_cutoff  = fget("Piscore_cutoff")

        self._validate()

    def _validate(self):
        errors = []
        s = self._stage

        if s == "1":
            if self.interpae_csv and not self.interpae_csv.exists():
                errors.append("Path_interpae_csv not found: {}".format(self.interpae_csv))

        if s == "2-4":
            for label, path in [
                ("Path_PDB_dir",        self.pdb_dir),
                ("Path_DeepTMHMM_dir",  self.deeptmhmm_dir),
                ("Path_DeepTMHMM_venv", self.deeptmhmm_py),
            ]:
                if path and not path.exists():
                    errors.append("{} not found: {}".format(label, path))

        if errors:
            raise ValueError(
                "Configuration errors:\n" +
                "\n".join("  - {}".format(e) for e in errors)
            )

        self.output_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _read_fasta(path):
        lines = path.read_text(encoding="utf-8").splitlines()
        seq = "".join(l.strip() for l in lines if not l.startswith(">"))
        if not seq:
            raise ValueError("FASTA file appears empty: {}".format(path))
        return seq

    def __repr__(self):
        lines = ["Config(stage={})".format(self._stage)]
        lines.append("  output_dir       = {}".format(self.output_dir))
        if self.interpae_csv:
            lines.append("  interpae_csv     = {}".format(self.interpae_csv))
        if self.pdb_dir:
            lines.append("  pdb_dir          = {}".format(self.pdb_dir))
        if self.deeptmhmm_dir:
            lines.append("  deeptmhmm_dir    = {}".format(self.deeptmhmm_dir))
        lines += [
            "  iptm_cutoff      = {}".format(self.iptm_cutoff),
            "  dockq_cutoff     = {}".format(self.dockq_cutoff),
            "  piscore_cutoff   = {}".format(self.piscore_cutoff),
            "  distance_cutoff  = {} A".format(self.distance_cutoff),
            "  pae_cutoff       = {}".format(self.pae_cutoff),
            "  max_workers_gpu  = {}".format(self.max_workers_gpu),
            "  max_workers_cpu  = {}".format(self.max_workers_cpu),
        ]
        return "\n".join(lines)
