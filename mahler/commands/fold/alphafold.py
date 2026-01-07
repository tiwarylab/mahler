"""AlphaFold orchestration helpers."""

from __future__ import annotations

from pathlib import Path
from os import PathLike
from mahler.utils.fasta import parse_fasta

__all__ = ["predict", "post_process"]

try:  # pragma: no cover - optional dependency
    import af2rave.alphafold as af2fold
except ModuleNotFoundError as exc:  # pragma: no cover
    af2fold = None  # type: ignore[assignment]
    _IMPORT_ERROR = exc
else:
    _IMPORT_ERROR = None


def predict(
        sequence: PathLike, 
        template: PathLike, 
        output_prefix: PathLike = Path("alphafold")
) -> None:
    """Run AlphaFold predictions across a preset panel of MSA depths."""

    # check modules
    if af2fold is None:
        raise ModuleNotFoundError(
            "af2rave is required for folding support. Install `af2rave` first."
        ) from _IMPORT_ERROR

    # load sequence
    if Path(sequence).is_file():
        records = parse_fasta(sequence)
        if len(records) != 1:
            raise ValueError("Input FASTA must contain exactly one sequence")
        sequence = next(iter(records.values()))
    else:
        # try treat sequence as raw
        sequence = str(sequence).strip()

    af = af2fold.AlphaFold.from_sequence(sequence, name="agab")
    params = {
        "num_seeds": 128,
        "num_recycles": 5,
        "use_templates": True,
        "custom_template_path": str(template),
        "save_recycles": True,
    }

    msas = {
        "8:16": "msa_8_16", 
        "16:32": "msa_16_32", 
        "32:64": "msa_32_64", 
        "520:5120": "msa_full"
    }

    for msa_range, dir_name in msas.items():
        output_path = Path(output_prefix) / dir_name
        output_path.mkdir(parents=True, exist_ok=True)
        af.predict(msa=msa_range, output_dir=str(output_path), **params)

def post_process(
        output_prefix: PathLike = Path("alphafold"), 
        destination: PathLike = Path("structures")
) -> None:
    """Collect AlphaFold structures into a single directory with readable names."""
    base = Path(output_prefix)
    if not base.exists():
        raise FileNotFoundError(f"{base} does not exist")

    structures_dir = base / destination
    structures_dir.mkdir(parents=True, exist_ok=True)

    for subdir in base.iterdir():
        if not subdir.is_dir() or not subdir.name.startswith("msa_"):
            continue

        prefix = subdir.name
        for pdb_file in subdir.glob("*.pdb"):
            name = pdb_file.name
            model = _extract_segment(name, "model_")
            seed = _extract_segment(name, "seed_", ".")
            recycle = _extract_segment(name, ".r", ".pdb", digits_only=True)
            if recycle is None:
                continue

            print(f"[{prefix}] model: {model}, seed: {seed}, recycle: {recycle}")
            destination = structures_dir / f"{prefix}_m.{model}_s.{seed}_r.{recycle}.pdb"
            destination.write_bytes(pdb_file.read_bytes())


def _extract_segment(
        filename: str, 
        start_marker: str, 
        end_marker: str = "_", 
        digits_only: bool = False
):
    start = filename.find(start_marker)
    if start == -1:
        return None
    start += len(start_marker)
    end = filename.find(end_marker, start)
    if end == -1:
        end = len(filename)
    value = filename[start:end]
    if digits_only and (not value.isdigit() or len(value) != 1):
        return None
    return value
