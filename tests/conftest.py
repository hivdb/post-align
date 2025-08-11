"""Pytest configuration for tests requiring third-party shims."""

from __future__ import annotations

import sys
import types
from collections.abc import Generator, Iterable, Iterator
from pathlib import Path
from unittest.mock import patch
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


@pytest.fixture(autouse=True)
def rich_shim() -> Generator[None, None, None]:
    """Provide a minimal :mod:`rich` shim with ``print``."""
    fake_rich = types.ModuleType("rich")
    fake_rich.print = print  # type: ignore[attr-defined]
    with patch.dict(sys.modules, {"rich": fake_rich}):
        yield


@pytest.fixture(autouse=True)
def more_itertools_shim() -> Generator[None, None, None]:
    """Provide a minimal :mod:`more_itertools` shim."""
    import itertools

    def chunked(iterable: Iterable[Any], n: int) -> Iterator[tuple[Any, ...]]:
        iterator = iter(iterable)
        while True:
            chunk = tuple(itertools.islice(iterator, n))
            if not chunk:
                break
            yield chunk

    fake_more = types.ModuleType("more_itertools")
    fake_more.chunked = chunked  # type: ignore[attr-defined]
    with patch.dict(sys.modules, {"more_itertools": fake_more}):
        yield


@pytest.fixture(autouse=True)
def orjson_shim() -> Generator[None, None, None]:
    """Provide a minimal :mod:`orjson` shim using :mod:`json`."""
    import json

    fake_orjson = types.ModuleType("orjson")
    fake_orjson.OPT_INDENT_2 = 2  # type: ignore[attr-defined]

    def dumps(obj: Any, *, option: int | None = None) -> bytes:
        indent = 2 if option == fake_orjson.OPT_INDENT_2 else None
        return json.dumps(obj, indent=indent).encode()

    fake_orjson.dumps = dumps  # type: ignore[attr-defined]
    with patch.dict(sys.modules, {"orjson": fake_orjson}):
        yield
