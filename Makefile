# Use the same Python as CI's lint job and Read the Docs, so local results match CI.
# Override with e.g. `make test UV_PYTHON=3.13`.
export UV_PYTHON ?= 3.12

# Every uv command selects the same extras as CI. Otherwise, when uv rebuilds .venv
# (e.g. after UV_PYTHON changes), the optional packages would be silently dropped.
UV_EXTRAS := --extra optional --extra property_prediction

.PHONY: install pre-commit test docs ci-local

install:
	uv sync $(UV_EXTRAS) --dev

pre-commit:
	uv run $(UV_EXTRAS) pre-commit run --all-files

test:
	uv run $(UV_EXTRAS) pytest --cov=smact -v

# Same settings as Read the Docs (.readthedocs.yaml: fail_on_warning, not nitpicky).
# -E rereads every source file, so warnings from an earlier run are never hidden.
docs:
	cd docs && uv run $(UV_EXTRAS) sphinx-build -E -W --keep-going -b html . _build/html/

ci-local: pre-commit test docs
