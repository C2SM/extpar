import re
import os
import sys
import argparse
from pathlib import Path

def find_markdown_files(base_path):
    """Recursively find all markdown files in the given directory."""
    return Path(base_path).rglob('*.md')

def extract_relative_links(content):
    """Extract relative links from markdown content."""
    # This regex matches markdown links. It captures the link target.
    pattern = re.compile(r'\[.*?\]\((?!http)(.*?\.md)\)')
    return pattern.findall(content)

def check_links(file_path, markdown_files):
    """Check if the relative links in a file are valid."""
    invalid_links = []
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()
        relative_links = extract_relative_links(content)
        for link in relative_links:
            # Construct the absolute path of the linked file
            linked_file = (file_path.parent / link).resolve()
            if not linked_file.is_file():
                invalid_links.append(link)
    return invalid_links

def main(base_path):
    markdown_files = set(find_markdown_files(base_path))
    all_invalid_links = {}

    for md_file in markdown_files:
        invalid_links = check_links(md_file)
        if invalid_links:
            all_invalid_links[md_file] = invalid_links

    if all_invalid_links:
        print("Found invalid links in the following files:")
        for file, links in all_invalid_links.items():
            print(f"{file}:")
            for link in links:
                print(f"  - {link}")
    else:
        print("No invalid links found.")

def main(base_path):
    markdown_files = set(find_markdown_files(base_path))
    all_invalid_links = {}
    all_checked_links = set()  # Set to keep track of all checked links

    for md_file in markdown_files:
        invalid_links = check_links(md_file, markdown_files)
        for link in extract_relative_links(open(md_file, 'r', encoding='utf-8').read()):
            all_checked_links.add((md_file, link))  # Add every checked link
        if invalid_links:
            all_invalid_links[md_file] = invalid_links

    print("Checked the following links:")
    for file, link in all_checked_links:
        print(f"{file}: {link}")

    if all_invalid_links:
        print("\nFound invalid links in the following files:")
        for file, links in all_invalid_links.items():
            print(f"{file}:")
            for link in links:
                print(f"  - {link}")
        sys.exit(1)  # Exit with a non-zero exit code if there are broken links
    else:
        print("\nNo invalid links found.")
        sys.exit(0)  # Exit with 0 if successful

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check markdown links in files.')
    parser.add_argument('-p', '--path', default=os.getcwd(), help='Base path to search for markdown files. Defaults to current working directory.')
    args = parser.parse_args()

    main(args.path)