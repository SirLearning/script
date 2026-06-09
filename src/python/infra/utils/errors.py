"""Shared error helpers for infra I/O and analytics modules."""

from __future__ import annotations

import logging
import sys

logger = logging.getLogger(__name__)


class DataLoadError(RuntimeError):
    """Raised when an input file cannot be read or parsed."""


def configure_logging(level: int = logging.INFO) -> None:
    """Configure root logging once for CLI / Nextflow-invoked scripts."""
    if logging.getLogger().handlers:
        return
    logging.basicConfig(
        level=level,
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )


def fail(message: str, *, exc: type[Exception] = DataLoadError) -> None:
    """Log an error and raise so callers (including Nextflow) get a non-zero exit."""
    logger.error(message)
    raise exc(message)
