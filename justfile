default: lint format check test

# Install dependencies
install:
    uv sync

# Run linting checks
lint:
    uv run ruff check src

# Format code
format:
	uv run ruff check --select I --fix src
	uv run ruff format src

# Run ty type checker
ty:
    uv run ty check src

# Run type checking
check:
    just lint
    just ty
    just test

# Run tests
test:
    uv run pytest tests

# Run tests with verbose output
test-v:
    uv run pytest tests -v

# Run tests for a specific file
test-file FILE:
    uv run pytest {{FILE}} -v

# Run tests with coverage (terminal)
cov:
    uv run pytest tests --cov=src/fastatacular --cov-report=term-missing

# Run tests with coverage (XML for Codecov) + JUnit XML for test results
test-cov:
    uv run pytest tests --cov=src/fastatacular --cov-report=xml --junitxml=junit.xml -o junit_family=legacy

# Remove cache and compiled files
clean:
    find . -type d -name __pycache__ -exec rm -rf {} +
    find . -type d -name .pytest_cache -exec rm -rf {} +
    find . -type d -name .ruff_cache -exec rm -rf {} +
    find . -name "*.pyc" -delete

# Build the package
build:
    uv build

# Install dev dependencies (alias for install)
dev:
    uv sync

# --- release (standard tacular-omics recipes; canonical copy in the workspace templates/) ---

# Set the version everywhere and date the [Unreleased] changelog section
set-version version:
    python scripts/release_version.py sync --set {{version}}

# Copy __version__ to CITATION.cff / .zenodo.json after editing it by hand
sync-version:
    python scripts/release_version.py sync

# Fail if version metadata disagrees
check-version:
    python scripts/release_version.py check
