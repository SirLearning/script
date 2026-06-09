"""Emit once-per-process deprecation warnings for genetics.wheat shims."""

from __future__ import annotations

import warnings

_WARNED: set[str] = set()


def warn_deprecated_shim(module_name: str, canonical: str) -> None:
    if module_name in _WARNED:
        return
    _WARNED.add(module_name)
    warnings.warn(
        f"genetics.wheat.{module_name} is deprecated; import from {canonical} instead.",
        DeprecationWarning,
        stacklevel=3,
    )
