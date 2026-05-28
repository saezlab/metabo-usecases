"""Unified logging helper for Python pipeline scripts.

Writes to the file pointed at by ``METABO_FIGURES_LOG`` using the line
format documented in
``specs/001-figures-pipeline/contracts/log-format.md``:

    <ISO-8601-with-tz> [Python][<level>][<component>] <message>

Falls back to ``logs/orphan-<PID>.log`` with a stderr warning when the
env var is missing.
"""

from __future__ import annotations

import logging
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


class _IsoTzFormatter(logging.Formatter):
    """Format times as ISO-8601 with the local timezone offset."""

    def formatTime(self, record: logging.LogRecord, datefmt: str | None = None) -> str:
        dt = datetime.fromtimestamp(record.created).astimezone()
        return dt.strftime("%Y-%m-%dT%H:%M:%S%z")


def _resolve_log_path() -> Path:
    raw = os.environ.get("METABO_FIGURES_LOG")
    if raw:
        return Path(raw)

    fallback = Path("logs") / f"orphan-{os.getpid()}.log"
    fallback.parent.mkdir(parents=True, exist_ok=True)
    os.environ["METABO_FIGURES_LOG"] = str(fallback)
    print(
        f"[python/utils/log.py] WARNING: METABO_FIGURES_LOG not set; "
        f"using {fallback}",
        file=sys.stderr,
    )
    return fallback


def configure_pipeline_log(component: str) -> logging.Logger:
    """Configure the root logger to append to the unified pipeline log.

    Parameters
    ----------
    component
        Logical owner string for this script's log lines. Conventionally
        one of ``svg_assemble:<id>``, ``cosmos_pkn``, ``check:<name>``,
        etc.

    Returns
    -------
    logging.Logger
        The configured component logger; use it for all subsequent log
        calls (``logger.info(...)``, ``logger.warning(...)``).
    """
    log_path = _resolve_log_path()

    handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    handler.setFormatter(
        _IsoTzFormatter(
            fmt=f"%(asctime)s [Python][%(levelname)s][{component}] %(message)s",
        )
    )

    root = logging.getLogger()
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(logging.INFO)

    return logging.getLogger(component)
