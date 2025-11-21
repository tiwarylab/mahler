"""Implementation of `mahler info`."""

from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from importlib import import_module

from .. import __version__
from .base import CommandResult


def _module_version(module_name: str, attr: str = "__version__") -> str:
    try:
        module = import_module(module_name)
    except ModuleNotFoundError:
        return "Not Found"

    value = getattr(module, attr, None)
    if value:
        return str(value)

    if module_name == "openmm":
        version_module = getattr(module, "version", None)
        if version_module:
            alt = getattr(version_module, "full_version", None) or getattr(
                version_module,
                "short_version",
                None,
            )
            if alt:
                return str(alt)
    return "Not Found"


def _torch_versions() -> tuple[str, str]:
    try:
        module = import_module("torch")
    except ModuleNotFoundError:
        return "Not Found", "Not Found"

    torch_version = getattr(module, "__version__", "Unknown")
    cuda_version = getattr(module.version, "cuda", None) if hasattr(module, "version") else None
    if callable(cuda_version):
        cuda_version = cuda_version()
    cuda_value = cuda_version or "Unavailable"
    return str(torch_version), str(cuda_value)


class InfoCommand:
    """Display basic diagnostic information about the CLI environment."""

    name = "info"
    help = "Display MAHLER CLI and environment information."

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--json",
            action="store_true",
            help="Print machine-readable JSON output.",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        torch_version, cuda_version = _torch_versions()
        payload = {
            "name": "mahler",
            "version": __version__,
            "platform": platform.platform(),
            "node": platform.node(),
            "python_version": platform.python_version(),
            "timestamp": datetime.now(tz=timezone.utc).isoformat(),
            "cuda_version": cuda_version,
            "pytorch_version": torch_version,
            "openmm_version": _module_version("openmm"),
            "af2rave_version": _module_version("af2rave"),
        }
        if args.json:
            print(json.dumps(payload, indent=2))
        else:
            print(f"MAHLER version {payload['version']}")
            print(f"Platform: {payload['platform']}")
            print(f"Node:     {payload['node']}")
            print(f"Python:   {payload['python_version']}")
            print(f"PyTorch:  {payload['pytorch_version']}")
            print(f"CUDA:     {payload['cuda_version']}")
            print(f"OpenMM:   {payload['openmm_version']}")
            print(f"af2rave:  {payload['af2rave_version']}")
        return CommandResult().exit_code
