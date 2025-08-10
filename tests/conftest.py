"""Pytest configuration for tests requiring third-party shims."""

from __future__ import annotations

import sys
import types
from collections.abc import Generator
from pathlib import Path
from unittest.mock import patch
from typing import Any, Callable

import pytest


@pytest.fixture(autouse=True)
def cython_shim() -> Generator[None, None, None]:
    """Provide a dummy `cython` module for tests.

    Some utilities import :mod:`cython` for type hints. The actual Cython
    package is not required for these tests, so a lightweight shim is
    injected into :data:`sys.modules` to satisfy the import.
    """
    fake_cython = types.ModuleType("cython")
    fake_cython.int = int  # type: ignore[attr-defined]
    fake_cython.ccall = lambda func: func  # type: ignore[attr-defined]
    fake_cython.cfunc = lambda func: func  # type: ignore[attr-defined]
    fake_cython.inline = lambda func: func  # type: ignore[attr-defined]
    fake_cython.cclass = lambda cls: cls  # type: ignore[attr-defined]
    fake_cython.bint = bool  # type: ignore[attr-defined]

    def declare(
        type_: object,
        value: object | None = None,
        **_: object,
    ) -> object:
        if value is not None:
            return value
        if type_ is int:
            return 0
        if type_ is bool:
            return False
        return None

    fake_cython.declare = declare  # type: ignore[attr-defined]

    def _returns(_: object) -> Callable[[Any], Any]:  # type: ignore[misc]
        return lambda func: func

    fake_cython.returns = _returns  # type: ignore[attr-defined]

    project_root = str(Path(__file__).resolve().parents[1])
    new_path = [project_root] + sys.path
    with patch.dict(sys.modules, {"cython": fake_cython}), patch.object(
        sys, "path", new_path
    ):
        yield


@pytest.fixture(autouse=True)
def pafpy_shim() -> Generator[None, None, None]:
    """Provide a minimal :mod:`pafpy` shim for parser imports."""
    fake_pafpy = types.ModuleType("pafpy")
    fake_pafpy.PafRecord = object  # type: ignore[attr-defined]
    fake_pafpy.Strand = object  # type: ignore[attr-defined]
    with patch.dict(sys.modules, {"pafpy": fake_pafpy}):
        yield
