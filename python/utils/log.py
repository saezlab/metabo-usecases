"""
Logging of Python steps into the unified pipeline log.

Lines are appended to the file named by ``METABO_FIGURES_LOG`` (set by
``rebuild.R``) in the format the R side writes (``R/utils-log.R``):

    <ISO-8601 time with offset> [Python][<level>][<component>] <message>

Without ``METABO_FIGURES_LOG``, the log goes to
``build/logs/orphan-<PID>.log`` in the repository, with a warning on stderr.
"""

from __future__ import annotations

import logging
import os
import sys
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


class _IsoTzFormatter(logging.Formatter):
    """
    Formats times as ISO-8601 with the local time zone offset.
    """

    def formatTime(
        self,
        record: logging.LogRecord,
        datefmt: str | None = None,
    ) -> str:

        dt = datetime.fromtimestamp(record.created).astimezone()

        return dt.strftime('%Y-%m-%dT%H:%M:%S%z')


def _log_path() -> Path:
    """
    The log file of this run, or a fallback in ``build/logs/``.
    """

    raw = os.environ.get('METABO_FIGURES_LOG')

    if raw:

        return Path(raw)

    fallback = REPO_ROOT / 'build' / 'logs' / f'orphan-{os.getpid()}.log'
    fallback.parent.mkdir(parents = True, exist_ok = True)
    os.environ['METABO_FIGURES_LOG'] = str(fallback)
    print(
        '[python/utils/log.py] WARNING: METABO_FIGURES_LOG not set; '
        f'using {fallback}',
        file = sys.stderr,
    )

    return fallback


def configure_pipeline_log(component: str) -> logging.Logger:
    """
    Sends the log records of this process to the pipeline log.

    Call it in the ``if __name__ == '__main__':`` block of a script, before
    anything else logs.

    Args:
        component: Owner of the log lines, e.g. ``analysis:<name>`` or
            ``check:<name>``.

    Returns:
        The logger of the component; use it for all further messages.
    """

    handler = logging.FileHandler(_log_path(), mode = 'a', encoding = 'utf-8')
    handler.setFormatter(
        _IsoTzFormatter(
            fmt = f'%(asctime)s [Python][%(levelname)s][{component}] '
            '%(message)s',
        ),
    )

    root = logging.getLogger()
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(logging.INFO)

    return logging.getLogger(component)
