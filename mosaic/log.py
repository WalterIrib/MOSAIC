"""
mosaic/log.py
Centralised logging for MOSAIC. Writes to both stderr (coloured if tty) and
a rotating session.log file in the project WD.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

# ─────────────────────────────────────────────────────────────────────────────
# ANSI colours (only when writing to a real terminal)
# ─────────────────────────────────────────────────────────────────────────────

_USE_COLOR = sys.stderr.isatty()

_COL = {
    "reset":  "\033[0m",
    "bold":   "\033[1m",
    "grey":   "\033[90m",
    "blue":   "\033[94m",
    "green":  "\033[92m",
    "yellow": "\033[93m",
    "red":    "\033[91m",
    "purple": "\033[95m",
}


def _col(text: str, *codes: str) -> str:
    if not _USE_COLOR:
        return text
    prefix = "".join(_COL[c] for c in codes if c in _COL)
    return f"{prefix}{text}{_COL['reset']}"


# ─────────────────────────────────────────────────────────────────────────────
# Custom formatter
# ─────────────────────────────────────────────────────────────────────────────

class MosaicFormatter(logging.Formatter):
    _ICONS = {
        logging.DEBUG:    ("DBG", "grey"),
        logging.INFO:     ("INF", "blue"),
        logging.WARNING:  ("WRN", "yellow"),
        logging.ERROR:    ("ERR", "red"),
        logging.CRITICAL: ("CRT", "red"),
    }

    def format(self, record: logging.LogRecord) -> str:
        icon, colour = self._ICONS.get(record.levelno, ("???", "grey"))
        tag = _col(f"[{icon}]", colour, "bold")
        msg = super().format(record)
        # Replace the default levelname prefix with our tag
        return f"{tag} {msg}"


# ─────────────────────────────────────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────────────────────────────────────

_logger: Optional[logging.Logger] = None


def setup_logging(log_dir: Optional[Path] = None,
                  level: int = logging.INFO) -> logging.Logger:
    """
    Initialise MOSAIC logging.
    Call once from the CLI entrypoint.
    """
    global _logger
    logger = logging.getLogger("mosaic")
    logger.setLevel(level)
    logger.handlers.clear()

    # stderr handler
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(MosaicFormatter("%(message)s"))
    logger.addHandler(sh)

    # file handler (when WD is known)
    if log_dir is not None:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_dir / "session.log", mode="a", encoding="utf-8")
        fh.setFormatter(logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        ))
        logger.addHandler(fh)

    _logger = logger
    return logger


def get_logger() -> logging.Logger:
    if _logger is None:
        return setup_logging()
    return _logger


# Convenience wrappers so callers don't import logging directly
def info(msg: str)     -> None: get_logger().info(msg)
def warn(msg: str)     -> None: get_logger().warning(msg)
def error(msg: str)    -> None: get_logger().error(msg)
def debug(msg: str)    -> None: get_logger().debug(msg)
def success(msg: str)  -> None: get_logger().info(_col(msg, "green", "bold"))
def section(msg: str)  -> None: get_logger().info(_col(f"\n{'─'*60}\n  {msg}\n{'─'*60}", "blue", "bold"))
