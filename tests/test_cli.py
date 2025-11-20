"""Smoke tests for the Tempo CLI backbone."""

from __future__ import annotations

import argparse

from tempo import cli
from tempo.commands import available_commands


def test_parser_registers_known_commands() -> None:
    parser = cli.build_parser(available_commands())
    assert isinstance(parser, argparse.ArgumentParser)

    args = parser.parse_args(["info"])
    assert getattr(args, "handler").name == "info"


def test_available_commands_listing() -> None:
    names = [command.name for command in available_commands()]
    assert names == ["info", "fold", "rave", "ifmetad", "reweight"]


def test_main_shows_help_without_command(capsys) -> None:
    exit_code = cli.main([])
    captured = capsys.readouterr()
    assert exit_code == 0
    assert "usage:" in captured.out.lower()
