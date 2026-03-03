.PHONY: install pre-commit test ci-local

install:
	uv sync --extra optional --extra property_prediction --dev

pre-commit:
	uv run pre-commit run --all-files

test:
	uv run pytest --cov=smact -v

ci-local: pre-commit test
