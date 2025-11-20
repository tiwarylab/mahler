"""Allow `python -m tempo` execution."""

from __future__ import annotations

from .cli import main


def run() -> int:
    """Delegate module execution to the CLI entry point."""
    return main()


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(run())
