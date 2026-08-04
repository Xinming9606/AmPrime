#!/usr/bin/env python3
"""Detect and install external AmPrime command-line dependencies."""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess

log = logging.getLogger(__name__)
REQUIRED_TOOLS = ("vsearch", "muscle")


def _log_process_output(completed: subprocess.CompletedProcess[str]) -> None:
    if completed.stdout:
        log.info(completed.stdout.rstrip())
    if completed.stderr:
        log.info(completed.stderr.rstrip())


def ensure_tool(tool: str) -> str:
    """Return a tool path, installing it with Scoop on Windows if needed."""
    executable = shutil.which(tool)
    if executable:
        return executable

    if os.name != "nt":
        raise RuntimeError(
            f"{tool} is required but was not found on PATH. Install it with "
            "Pixi/Conda and retry."
        )

    scoop = shutil.which("scoop")
    if not scoop:
        raise RuntimeError(
            f"{tool} is required on Windows, but Scoop was not found. "
            f"Install Scoop, then run 'scoop install {tool}'."
        )

    log.info("%s not found; installing it with Scoop", tool)
    completed = subprocess.run(  # noqa: S603 - executable came from PATH lookup.
        [scoop, "install", tool],
        capture_output=True,
        text=True,
        check=False,
    )
    _log_process_output(completed)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Scoop failed to install {tool} with exit code "
            f"{completed.returncode}. Run 'scoop install {tool}' manually."
        )

    executable = shutil.which(tool)
    if not executable:
        raise RuntimeError(
            f"Scoop reported a successful {tool} installation, but {tool} was "
            "not found on PATH. Restart the shell and retry."
        )
    return executable


def ensure_required_tools() -> dict[str, str]:
    """Ensure every external executable required by the workflow is present."""
    return {tool: ensure_tool(tool) for tool in REQUIRED_TOOLS}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    tools = ensure_required_tools()
    for name, executable in tools.items():
        log.info("%s: %s", name, executable)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
