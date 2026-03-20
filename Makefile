.DEFAULT_GOAL := help

help:  ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install-dev:  ## Install dev dependencies via pipenv
	pipenv install --dev

test:  ## Run test suite
	pipenv run python -m pytest tests/ -v

lint:  ## Run ruff check + mypy
	pipenv run ruff check postalign/ tests/
	pipenv run mypy postalign/

format:  ## Run ruff format
	pipenv run ruff format postalign/ tests/
	pipenv run ruff check --fix postalign/ tests/

benchmark:  ## Run performance benchmarks (saves JSON + generates charts)
	@mkdir -p benchmarks/data/charts
	pipenv run python -m pytest benchmarks/bench_codon_align.py -v \
		--benchmark-only \
		--benchmark-json=benchmarks/data/results.json
	pipenv run python benchmarks/generate_report.py

benchmark-report:  ## Re-generate charts from existing benchmark JSON
	pipenv run python benchmarks/generate_report.py

build-ext:  ## Build Cython extensions
	pipenv run python setup.py build_ext --inplace

clean:  ## Remove build artifacts, .so, __pycache__, etc.
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name '*.so' -delete 2>/dev/null || true
	find . -type f -name '*.c' -not -name 'codon_alignment.c' -delete 2>/dev/null || true
	rm -rf build/ dist/ *.egg-info .pytest_cache .mypy_cache

build-docker-builder:
	@docker pull ubuntu:18.04
	# @docker build . --no-cache -t hivdb/post-align-builder:latest
	@docker build . -t hivdb/post-align-builder:latest

push-docker-builder: build-docker-builder
	@docker push hivdb/post-align-builder:latest

inspect-builder:
	@docker run \
		--mount type=bind,source=$(PWD),target=/app/post-align \
		--workdir /app --rm -it \
		hivdb/post-align-builder:latest \
		bash

dist/linux-amd64: $(shell find postalign -type f -path "*.py" | sed 's#\([| ]\)#\\\1#g')
	@docker run \
		--mount type=bind,source=$(PWD),target=/app/post-align \
		--workdir /app --rm -it \
		hivdb/post-align-builder:latest \
		/build.sh

dist/postalign_linux-amd64.tar.gz: dist/linux-amd64 
	@cd dist/linux-amd64 && \
		rm -f mpostalign_linux-amd64.tar.gz && \
		tar zcf mpostalign_linux-amd64.tar.gz mpostalign && \
		mv mpostalign_linux-amd64.tar.gz ..

dist: dist/postalign_linux-amd64.tar.gz
build: dist

.PHONY: help install-dev test lint format benchmark benchmark-report build-ext clean \
	build-docker-builder push-docker-builder
