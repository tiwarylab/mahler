"""Simple helpers for working with FASTA files."""

from __future__ import annotations

from pathlib import Path
from typing import IO

__all__ = ["parse_fasta"]


def parse_fasta(source: str | Path | IO[str]) -> dict[str, str]:
    """Return a mapping of FASTA identifiers to their sequences.

    Parameters
    ----------
    source:
        Either a filesystem path or an already opened text handle. Empty lines
        are ignored and sequences are stripped of surrounding whitespace. Raises
        ``ValueError`` when the input is malformed (e.g. duplicate identifiers
        or sequence data before the first header).
    """

    def _flush(current_id: str | None, current_seq: list[str]) -> None:
        if current_id is None:
            return
        sequences[current_id] = "".join(current_seq)

    if hasattr(source, "read"):
        handle = source  # type: ignore[assignment]
        should_close = False
    else:
        handle = Path(source).expanduser().open("r", encoding="utf-8")
        should_close = True

    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []

    try:
        for lineno, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header on line {lineno}")
                if header in sequences:
                    raise ValueError(
                        f"Duplicate FASTA identifier '{header}' on line {lineno}"
                    )
                _flush(current_id, current_seq)
                current_id, current_seq = header, []
                continue

            if current_id is None:
                raise ValueError(
                    f"Encountered sequence data before first header on line {lineno}"
                )
            current_seq.append(line)

        _flush(current_id, current_seq)
    finally:
        if should_close:
            handle.close()

    return sequences
