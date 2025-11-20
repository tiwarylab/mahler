"""Reusable command registry for the TEMPO CLI."""

from __future__ import annotations

from typing import List, Sequence, Type

from .base import Command
from .fold import FoldCommand
from .ifmetad import IfMetadCommand
from .info import InfoCommand
from .rave import RaveCommand
from .reweight import ReweightCommand


def available_commands() -> Sequence[Type[Command]]:
    """Return all command classes to register with argparse."""
    commands: List[Type[Command]] = [
        InfoCommand,
        FoldCommand,
        RaveCommand,
        IfMetadCommand,
        ReweightCommand,
    ]
    return commands


__all__ = ["Command", "available_commands"]
