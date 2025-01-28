# Contributing

This is a quick guide on how to follow best practice and contribute smoothly to `SMACT`.

## Workflow

We follow the [GitHub flow](https://docs.github.com/en/get-started/using-github/github-flow), using
branches for new work and pull requests for verifying the work.

The steps for a new piece of work can be summarised as follows:

1. Push up or create [an issue](https://github.com/WMD-group/SMACT/issues).
2. Fork the repository, i.e. go to the [`SMACT` repository on GitHub](https://github.com/WMD-group/SMACT) and click the "Fork" button to create a copy of the repository under your own account.
3. Clone your fork of the repository to your local machine.

   ```bash
   git clone https://github.com/<your-username>/SMACT.git
   cd SMACT
   ```

4. Install the [uv package manager](https://github.com/astral-sh/uv).

   ```bash
   pip install uv # installs uv package manager
   ```

5. Create a virtual environment for SMACT.

   Using `uv`:

   ```bash
   uv create venv # A virtual environment will be created in the current directory
   source .venv/bin/activate # Activate the virtual environment
   ```

   N.B. A virtual environment can also be created using [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) for those more familiar with mamba/conda.

6. Install the package in editable mode with the development, documentation and optional dependencies.

   ```bash
   uv pip install -e ".[dev,docs,optional]"
   pre-commit install # Install pre-commit hooks
   ```

7. Create a branch from master, with a sensible name that relates to the issue.

   ```bash
   git checkout -b <branch-name> # should be run from the master branch
   ```

8. Do the work and commit changes to the branch. Push the branch
   regularly to GitHub to make sure no work is accidentally lost.

   ```bash
   git add <files>
   git commit -m "A message describing the changes"
   git push origin <branch-name>
   ```

9. Write or update unit tests for the code you work on.
10. When you are finished with the work, ensure that all of the unit tests pass on your own machine by running `pytest -v` in the root directory of the repository.

11. Open a pull request [on the pull request page](https://github.com/WMD-group/SMACT/pulls).
12. If nobody acknowledges your pull request promptly, feel free to poke one of the main developers into action.
13. For keeping your repository up to date with the master repository, you can add it as a remote to your local repository.

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

When developing locally, it is recommended to install the python packages for development and documentation.

```bash
pip install -e ".[dev,docs]" # Should be ran from the root of the repository
```

This will allow you to run the tests locally with pytest as described in the main README,
as well as run pre-commit hooks to automatically format python files with [ruff](https://docs.astral.sh/ruff/).
To install the pre-commit hooks (only needs to be done once):

```bash
pre-commit install
pre-commit run --all-files # optionally run hooks on all files
```

Pre-commit hooks will check all files when you commit changes, automatically fixing any files which are not formatted correctly. Those files will need to be staged again before re-attempting the commit.

To render the documentation locally, you can use the following command:

```bash
cd docs
sphinx-build -nW --keep-going -b html . _build/html/
```

You should then be able to view the documentation by opening `_build/html/index.html` in a web browser. You should inspect the documentation to ensure that it is correctly formatted, any docstrings are present and that the code examples are correct.
