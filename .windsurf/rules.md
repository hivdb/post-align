# Post-Align Project Rules

## Domain
HIV bioinformatics post-alignment toolkit. Refines pairwise/multiple alignment sequences with codon-aware gap placement using BLOSUM62 and IUPAC scoring.

## Code Style
- **Strict mypy**: config in `pyproject.toml` `[tool.mypy]`. All functions must have type annotations.
- **Ruff**: linting and formatting via `[tool.ruff]` in `pyproject.toml`. Run `ruff check` and `ruff format`.
- **Cython decorators**: hot-path functions use `@cython.cfunc`, `@cython.inline`, `@cython.returns()` for optional Cython compilation. These are no-ops in pure Python mode.
- **Python >= 3.13** required.
- **Type hints**: use lowercase builtins (`list`, `dict`, `set`, `tuple`) — do NOT import from `typing`. Use `X | Y` union syntax instead of `Union[X, Y]`, `X | None` instead of `Optional[X]`.
- **pipenv** is the mandatory dependency manager. Use `pipenv install --dev` and `pipenv run` for all commands.
- **Dependency management**: runtime dependencies live in `pyproject.toml` `[project.dependencies]`. `Pipfile` references the project via `post-align = {editable = true, path = "."}` so there's a single source of truth. Dev-only dependencies (pytest, ruff, mypy, etc.) stay in `Pipfile` `[dev-packages]`. No `requirements.txt`.

## Architecture
- **NAPosition**: core per-nucleotide data model (`postalign/models/na_position.py`). Has `notation` (int), `pos` (int), `flag` (PositionFlag), `is_gap` (bool), `payload`.
- **Sequence**: generic sequence container wrapping `List[NAPosition]` with modifier history.
- **RefSeqPair**: `Tuple[Sequence, Sequence]` — the standard ref+seq pair passed through the pipeline.
- **Processor pipeline**: CLI commands are chained processors (`postalign/processor.py`). Intermediate processors yield `Iterable[RefSeqPair]`, output processors yield `Iterable[str]`.
- **Codon alignment**: `postalign/processors/codon_alignment.py` — the primary algorithm. Public API is `codon_align()` and `parse_gap_placement_score()`.

## Testing
- **pytest + hypothesis** for all tests.
- Tests target the **public API** (`codon_align`, `parse_gap_placement_score`), not internal functions.
- All Phase 2 optimization tiers (T1–T4) must pass the same test suite.
- Run tests with `make test`.

## Maintenance Scripts
- **All maintenance scripts must have entrypoints from `Makefile`.**
- Available targets: `make test`, `make lint`, `make format`, `make benchmark`, `make benchmark-report`, `make build-ext`, `make install-dev`, `make clean`.
- All `make` targets use `pipenv run` internally.

## Performance Tiers
- **T1**: Original code (baseline)
- **T2**: Original code + algorithm improvements (Option D)
- **T3**: Cython byte-array rewrite + Option D
- **T4**: Rust/PyO3 core + Option D
- All tiers live as separate files in `postalign/processors/` with the same `codon_align()` API.
