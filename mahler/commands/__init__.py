"""Reusable command registry for the TEMPO CLI."""

from __future__ import annotations

from .base import Command
from .fold import FoldCommand
from .ifmetad import IfMetadCommand
from .info import InfoCommand
from .rave import RaveCommand
from .reweight import ReweightCommand
from .mdprep import MDPrepCommand
from .mdrun import MDRunCommand


def available_commands() -> list[Command]:
    """Return all command classes to register with argparse."""
    commands: list[Command] = [
        InfoCommand,
        FoldCommand,
        RaveCommand,
        IfMetadCommand,
        ReweightCommand,
        MDPrepCommand,
        MDRunCommand,
    ]
    return commands

__all__ = ["Command", "available_commands"]
