import re
import argparse

def remove_github_actions_entries(changelog_path):
    with open(changelog_path, 'r') as file:
        lines = file.readlines()

    with open(changelog_path, 'w') as file:
        skip = False
        for line in lines:
            if re.search(r'\[github-actions\[bot\]\]', line):
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