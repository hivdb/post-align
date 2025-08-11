"""Pytest configuration for tests requiring third-party shims."""

from __future__ import annotations

import sys
import types
from collections.abc import Generator, Iterable, Iterator
from pathlib import Path
from unittest.mock import MagicMock, patch
from typing import Any, Callable, TextIO

import pytest

fake_typer = types.ModuleType("typer")


class _Typer:
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        pass

    def result_callback(
        self, func: Callable[..., Any] | None = None
    ) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        def decorator(inner: Callable[..., Any]) -> Callable[..., Any]:
            return inner

        if func is not None:
            return decorator(func)
        return decorator

    def callback(
        self, *args: Any, **kwargs: Any
    ) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            return func

        return decorator

    def command(
        self, *args: Any, **kwargs: Any
    ) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            return func

        return decorator


class _BadParameter(Exception):
    pass


fake_typer.Typer = _Typer  # type: ignore[attr-defined]
fake_typer.BadParameter = _BadParameter  # type: ignore[attr-defined]
fake_typer.FileText = TextIO  # type: ignore[attr-defined]
fake_typer.FileTextWrite = TextIO  # type: ignore[attr-defined]


class _Context:
    def __init__(self) -> None:
        self.params: dict[str, Any] = {}


class _CallbackParam:
    def __init__(self, name: str | None = None) -> None:
        self.name = name


def _option(*_args: Any, **_kwargs: Any) -> None:  # noqa: D401
    """Return a no-op Typer option placeholder."""
    return None


fake_typer.Context = _Context  # type: ignore[attr-defined]
fake_typer.CallbackParam = _CallbackParam  # type: ignore[attr-defined]
fake_typer.Option = _option  # type: ignore[attr-defined]
fake_typer.Argument = _option  # type: ignore[attr-defined]
sys.modules.setdefault("typer", fake_typer)


fake_cython = types.ModuleType("cython")
fake_cython.int = int  # type: ignore[attr-defined]
fake_cython.ccall = lambda func: func  # type: ignore[attr-defined]
fake_cython.cfunc = lambda func: func  # type: ignore[attr-defined]
fake_cython.inline = lambda func: func  # type: ignore[attr-defined]
fake_cython.cclass = lambda cls: cls  # type: ignore[attr-defined]
fake_cython.bint = bool  # type: ignore[attr-defined]


def _declare(type_: object, value: object | None = None, **_: object) -> object:
    """Mimic :func:`cython.declare` for test shimming."""
    if value is not None:
        return value
    if type_ is int:
        return 0
    if type_ is bool:
        return False
    return None


fake_cython.declare = _declare  # type: ignore[attr-defined]


def _returns(_: object) -> Callable[[Any], Any]:  # type: ignore[misc]
    return lambda func: func


fake_cython.returns = _returns  # type: ignore[attr-defined]
sys.modules.setdefault("cython", fake_cython)


@pytest.fixture(autouse=True)
def cython_shim() -> Generator[None, None, None]:
    """Ensure the project root is importable during tests."""
    project_root = str(Path(__file__).resolve().parents[1])
    new_path = [project_root] + sys.path
    with patch.object(sys, "path", new_path):
        yield


fake_pafpy = types.ModuleType("pafpy")


class _PafRecord(types.SimpleNamespace):
    """Minimal stand-in for :class:`pafpy.PafRecord`."""

    @classmethod
    def from_str(cls, line: str) -> "_PafRecord":
        fields = line.rstrip().split("\t")
        tags = {
            tag.split(":", 2)[0]: types.SimpleNamespace(value=tag.split(":", 2)[2])
            for tag in fields[12:]
        }
        return cls(
            qname=fields[0],
            qlen=int(fields[1]),
            qstart=int(fields[2]),
            qend=int(fields[3]),
            strand=fields[4],
            tname=fields[5],
            tlen=int(fields[6]),
            tstart=int(fields[7]),
            tend=int(fields[8]),
            mlen=int(fields[9]),
            blen=int(fields[10]),
            mapq=int(fields[11]),
            tags=tags,
        )


fake_pafpy.PafRecord = _PafRecord  # type: ignore[attr-defined]
fake_pafpy.Strand = types.SimpleNamespace(Reverse="-")  # type: ignore[attr-defined]
sys.modules.setdefault("pafpy", fake_pafpy)


@pytest.fixture(autouse=True)
def pafpy_shim() -> Generator[None, None, None]:
    """No-op fixture retaining the :mod:`pafpy` shim."""
    yield


fake_rich = types.ModuleType("rich")
fake_rich.print = MagicMock(side_effect=print)  # type: ignore[attr-defined]
sys.modules.setdefault("rich", fake_rich)


@pytest.fixture(autouse=True)
def rich_shim() -> Generator[None, None, None]:
    """Retain the :mod:`rich` shim for tests."""
    yield


import itertools


def _chunked(iterable: Iterable[Any], n: int) -> Iterator[tuple[Any, ...]]:
    """Simplified replacement for :func:`more_itertools.chunked`."""
    iterator = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(iterator, n))
        if not chunk:
            break
        yield chunk


fake_more = types.ModuleType("more_itertools")
fake_more.chunked = _chunked  # type: ignore[attr-defined]
sys.modules.setdefault("more_itertools", fake_more)


@pytest.fixture(autouse=True)
def more_itertools_shim() -> Generator[None, None, None]:
    """No-op fixture retaining the :mod:`more_itertools` shim."""
    yield


import json


fake_orjson = types.ModuleType("orjson")
fake_orjson.OPT_INDENT_2 = 2  # type: ignore[attr-defined]


def _orjson_dumps(obj: Any, *, option: int | None = None) -> bytes:
    """Serialize ``obj`` using :mod:`json` as a stand-in for :func:`orjson.dumps`."""
    indent = 2 if option == fake_orjson.OPT_INDENT_2 else None
    return json.dumps(obj, indent=indent).encode()


fake_orjson.dumps = _orjson_dumps  # type: ignore[attr-defined]
sys.modules.setdefault("orjson", fake_orjson)


@pytest.fixture(autouse=True)
def orjson_shim() -> Generator[None, None, None]:
    """No-op fixture to keep the :mod:`orjson` shim in place."""
    yield
