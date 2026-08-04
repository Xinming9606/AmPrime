#!/usr/bin/env python3
"""Small shared helpers for AmPrime workflow command-line tools."""

from __future__ import annotations

import logging
from pathlib import Path

IUPAC_COMPLEMENT_TABLE = str.maketrans(
    "ACGTRYMKSWHBVDNacgtrymkswhbvdn", "TGCAYRKMSWDVBHNtgcayrkmswdvbhn"
)


def configure_logging(log_path: str) -> None:
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def config_param(cli_value, cfg: dict, key: str):
    return cli_value if cli_value is not None else cfg.get(key)


def required_param(name: str, value):
    if value is None:
        raise SystemExit(
            f"missing --{name.replace('_', '-')} or config setting: {name}"
        )
    return value


def reverse_complement(seq: str) -> str:
    return seq.translate(IUPAC_COMPLEMENT_TABLE)[::-1]
