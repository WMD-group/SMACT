.PHONY: install lint typecheck test format ci-local

install:
	uv sync --extra optional --dev

lint:
	uv run ruff check smact/
	uv run ruff format --check smact/
	uv run codespell

typecheck:
	uv run pyright smact/

test:
	uv run pytest --cov=smact -v

format:
	uv run ruff check --fix smact/
	uv run ruff format smact/

ci-local: lint typecheck test
