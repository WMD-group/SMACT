"""Remove GitHub Actions bot entries from changelog."""

from __future__ import annotations

import argparse
import re


def remove_github_actions_entries(changelog_path: str) -> None:
    """Remove GitHub Actions bot entries from the changelog.

    Args:
        changelog_path: Path to the changelog file.
    """
    with open(changelog_path) as file:
        lines = file.readlines()

    with open(changelog_path, "w") as file:
        skip = False
        for line in lines:
            if re.search(r"\[github-actions\[bot\]\]", line):
                skip = True
            elif skip and line.strip() == "":
                skip = False
            elif not skip:
                file.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove GitHub Actions bot entries from changelog.")
    parser.add_argument("changelog_path", type=str, help="Path to the changelog file")
    args = parser.parse_args()

    remove_github_actions_entries(args.changelog_path)
