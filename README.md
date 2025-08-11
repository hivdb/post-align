# post-align

`post-align` is a Python toolkit and command-line interface for refining
pairwise or multiple sequence alignments. It provides utilities for parsing
common alignment formats and applying a pipeline of processors to generate
cleaned-up results and reports.

## Features

* Parsers for MSA, PAF and minimap2 outputs, including built-in execution of
  minimap2 to align raw reads against a reference.
* Codon-aware alignment and gap adjustment helpers.
* Processor pipeline with built-in steps such as `codon-alignment`,
  `trim-by-ref`, `save-fasta` and `save-json`.
* CLI entry point (`post-align`) for composing processors on the command line.

## Installation

```bash
pip install -e .
```

For development work, create the pipenv environment with

```bash
pipenv install --dev
```

## Usage

Run the command-line interface with an input file and a sequence of processors.
The typical workflow runs minimap2 to generate alignments and writes the
post‑processed result as JSON:

```bash
pipenv run post-align -i reads.fasta -r ref.fasta -o result.json -f MINIMAP2 \
    save-json
```

## Development

Linting, type checking and tests:

```bash
pipenv run flake8 .
pipenv run mypy .
pipenv run pytest tests/unit --cov=postalign --cov-report=term-missing
pipenv run behave tests/component  # downloads minimap2 2.17 automatically
```

The project currently relies on a minimal `setup.py` for building Cython
extensions while `pyproject.toml` provides the build backend. A future goal is
to consolidate all metadata into `pyproject.toml` once Cython support is
sufficient. To keep `requirements.txt` in sync with the Pipfile, regenerate it
after updating dependencies:

```bash
make requirements.txt
```

## License

This project is licensed under the [MIT License](LICENSE).
