# Use the same Python as CI's lint job and Read the Docs, so local results match CI.
# Override with e.g. `make test UV_PYTHON=3.13`.
export UV_PYTHON ?= 3.12

.PHONY: install pre-commit test docs ci-local

install:
	uv sync --extra optional --extra property_prediction --dev

pre-commit:
	uv run pre-commit run --all-files

test:
	uv run pytest --cov=smact -v

# Same settings as Read the Docs (.readthedocs.yaml: fail_on_warning, not nitpicky).
docs:
	cd docs && uv run sphinx-build -W --keep-going -b html . _build/html/

ci-local: pre-commit test docs
