"""Command protocol for TEMPO CLI subcommands."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Protocol


class Command(Protocol):
    """Protocol that each CLI command implementation must follow."""

    name: str
    help: str

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        """Add command-specific arguments."""

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        """Execute command logic."""


@dataclass
class CommandResult:
    """Simple value object for returning command execution details."""

    exit_code: int = 0
    message: str | None = None

    def emit(self) -> int:
        if self.message:
            print(self.message)
        return self.exit_code
