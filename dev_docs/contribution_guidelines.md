# Contributing to SMACT

This is a quick guide on how to follow best practice and contribute smoothly to `SMACT`.

## Workflow

We follow the [GitHub flow]
(<https://guides.github.com/introduction/flow/index.html>), using
branches for new work and pull requests for verifying the work.

The steps for a new piece of work can be summarised as follows:

1. Push up or create [an issue](https://guides.github.com/features/issues).
2. Create a branch from master, with a sensible name that relates to the issue.
3. Do the work and commit changes to the branch. Push the branch
   regularly to GitHub to make sure no work is accidentally lost.
4. Write or update unit tests for the code you work on.
5. When you are finished with the work, ensure that all of the unit
   tests pass on your own machine.
6. Open a pull request [on the SMACT pull request page](https://github.com/WMD-group/SMACT/pulls).
7. If nobody acknowledges your pull request promptly, feel free to poke one of the main developers into action.

## Pull requests

For a general overview of using pull requests on GitHub look [in the GitHub docs](https://help.github.com/en/articles/about-pull-requests).

When creating a pull request you should:

* Ensure that the title succinctly describes the changes so it is easy to read on the overview page
* Reference the issue which the pull request is closing

Recommended reading: [How to Write the Perfect Pull Request](https://github.blog/2015-01-21-how-to-write-the-perfect-pull-request/)

## Dev requirements

When developing `SMACT` locally, it is recommended to install the python packages in `requirements-dev.txt`.

```bash
pip install -r requirements-dev.txt
```

This will allow you to run the tests locally with pytest as described in the main README,
as well as run pre-commit hooks to automatically format python files with isort and black.
To install the pre-commit hooks (only needs to be done once):

```bash
pre-commit install
pre-commit run --all-files # optionally run hooks on all files
```

Pre-commit hooks will check all files when you commit changes, automatically fixing any files which are not formatted correctly. Those files will need to be staged again before re-attempting the commit.
