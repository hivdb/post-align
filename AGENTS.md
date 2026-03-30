# Agent Rules for post-align

## Code Quality

- **Strict ruff**: all code must pass `ruff check` and `ruff format` using the config in `pyproject.toml` `[tool.ruff]`. Target Python 3.13, line length 120, single-quote style.
- **Strict mypy**: all code must pass `mypy` using the config in `pyproject.toml` `[tool.mypy]`. Every function must have complete type annotations. Use lowercase builtins (`list`, `dict`, `set`, `tuple`), `X | Y` unions, and `X | None` optionals — never import from `typing`.
- **No `Any` or `object`** in type hints — always use explicit, concrete types. If a Rust/PyO3 extension type is involved, maintain a `.pyi` stub in `stubs/` so mypy resolves it properly.
- Run `make lint` to verify both ruff and mypy before considering any change complete.

## Common Commands

All commands use `pipenv run` internally. Run them via `make`:

| Command                | Description                                      |
|------------------------|--------------------------------------------------|
| `make install-dev`     | Install dev dependencies via pipenv              |
| `make test`            | Run test suite (`pytest tests/ -v`)              |
| `make lint`            | Run `ruff check` + `mypy`                        |
| `make format`          | Run `ruff format` + `ruff check --fix`           |
| `make build-ext`       | Build Cython extensions                          |
| `make benchmark`       | Run benchmarks, save JSON, generate charts       |
| `make benchmark-report`| Re-generate charts from existing benchmark JSON  |
| `make clean`           | Remove build artifacts, `.so`, `__pycache__`, etc|

## Testing

- Framework: **pytest** (with **hypothesis** for property-based tests).
- Test directory: `tests/` — mirrors the `postalign/` package structure.
- Test files are prefixed with `test_`. Examples:
  - `postalign/processors/codon_alignment.py` → `tests/processors/test_codon_alignment.py`
  - `postalign/models/na_position.py` → `tests/models/test_na_position.py`
  - `postalign/utils/cigar.py` → `tests/utils/test_cigar.py`
- Shared fixtures go in `tests/conftest.py`.
- Run tests with `make test`.

## Temporary Scripts

**Never feed large chunks of code to `python` or `bash` via stdin or `-c`.** The input gets truncated and the command never terminates. Instead, write temporary scripts to a file under `/tmp` and execute them:

```sh
# Good
pipenv run python /tmp/my_smoke_test.py

# Bad — will hang
pipenv run python -c "$(cat <<'EOF'
...large code block...
EOF
)"
```
