# Contributing

This is a quick guide on how to follow best practice and contribute smoothly to `SMACT`.

> **Note (v4.0.0):** As of v4.0.0, SMACT has adopted **GitHub Flow** as its
> branching strategy. The `develop` branch is no longer used. All feature
> branches should be created from `master`, and all pull requests should
> target `master`. Release versions are cut by tagging commits on `master`.

## Workflow

We follow [GitHub Flow](https://docs.github.com/en/get-started/using-github/github-flow), using
branches for new work and pull requests for verifying the work.

The steps for a new piece of work can be summarised as follows:

1. Push up or create [an issue](https://github.com/WMD-group/SMACT/issues).
2. Fork the repository, i.e. go to the [`SMACT` repository on GitHub](https://github.com/WMD-group/SMACT) and click the "Fork" button to create a copy of the repository under your own account.
3. Clone your fork of the repository to your local machine.

   ```bash
   git clone https://github.com/<your-username>/SMACT.git
   cd SMACT
   ```

4. Install the [uv package manager](https://docs.astral.sh/uv/getting-started/installation/).

   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

5. Install SMACT with all development and optional dependencies. This creates a virtual environment automatically and installs the package in editable mode.

   ```bash
   uv sync --extra optional --dev
   pre-commit install # Install pre-commit hooks
   ```

6. Create a branch from master, with a sensible name that relates to the issue.

   ```bash
   git checkout -b <branch-name> # should be run from the master branch
   ```

7. Do the work and commit changes to the branch. Push the branch
   regularly to GitHub to make sure no work is accidentally lost.

   ```bash
   git add <files>
   git commit -m "A message describing the changes"
   git push origin <branch-name>
   ```

8. Write or update unit tests for the code you work on.
9. When you are finished with the work, run the full CI pipeline locally to catch issues before pushing:

   ```bash
   make ci-local # runs pre-commit hooks and tests
   ```

10. Open a pull request [on the pull request page](https://github.com/WMD-group/SMACT/pulls).
11. If nobody acknowledges your pull request promptly, feel free to poke one of the main developers into action.
12. For keeping your repository up to date with the master repository, you can add it as a remote to your local repository.

```bash
git remote add upstream https://github.com/WMD-group/SMACT.git
```

Fetch the latest changes from the master branch to keep it up to date (make sure you are on the master branch ).

```bash
git checkout master
git pull upstream master
```

Reminder, `pull` is a combination of `fetch` and `merge`. This could lead to merge conflicts if you have made changes to the master branch.

## Pull requests

For a general overview of using pull requests on GitHub look [in the GitHub docs](https://help.github.com/en/articles/about-pull-requests).

When creating a pull request you should:

- Ensure that the title succinctly describes the changes so it is easy to read on the overview page
- Reference the issue which the pull request is closing

Recommended reading: [How to Write the Perfect Pull Request](https://github.blog/2015-01-21-how-to-write-the-perfect-pull-request/)

## Dev requirements

All development dependencies are managed in `pyproject.toml`. To install them:

```bash
uv sync --extra optional --dev
```

### Pre-commit hooks

Pre-commit hooks automatically run [ruff](https://docs.astral.sh/ruff/) (lint and format), [codespell](https://github.com/codespell-project/codespell), and [pyright](https://github.com/microsoft/pyright) when you commit changes. To install (only needs to be done once):

```bash
pre-commit install
pre-commit run --all-files # optionally run hooks on all files
```

Files auto-fixed by hooks will need to be re-staged before re-attempting the commit.

### Makefile targets

The Makefile provides convenient shortcuts for common tasks:

```bash
make install    # install all dependencies
make pre-commit # run all pre-commit hooks (ruff, pyright, codespell, prettier, etc.)
make test       # run pytest with coverage
make ci-local   # run pre-commit + test (matches CI)
```

### Documentation

To render the documentation locally:

```bash
cd docs
sphinx-build -nW --keep-going -b html . _build/html/
```

Then open `_build/html/index.html` in a browser to inspect formatting, docstrings, and code examples.
